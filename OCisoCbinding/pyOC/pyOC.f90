! This module is simply a 'sanitized' and f90wrapp-compatible version of the different subroutines available in the OCASI interface (liboctq.F90)
! all subroutines of liboctq.F90 'wrapped' into a new subroutine defined here (the subroutine name is the same with a 'py' prefix)
!
module RawOpenCalphad
  use liboctq
  implicit none
  
  type eq_wrapper
    type(gtp_equilibrium_data), pointer :: ceq
  end type eq_wrapper
  
  type comp_wrapper
    integer :: n
    character(24) :: compnames(maxc)
  end type comp_wrapper
  
!  integer, parameter :: maxel=100,maxsp=1000,maxph=600,maxsubl=10,maxconst=1000
  
contains

  function pygeterr() result(errorcode)
    integer :: errorcode
    errorcode = gx%bmperr
  end function pygeterr
  
  subroutine pyseterr(errorcode) 
    integer, intent(in) :: errorcode
    gx%bmperr = errorcode
  end subroutine pyseterr

  subroutine pytqini(n,eq)
    implicit none
    integer :: n
    type(eq_wrapper), intent(out) :: eq
    call tqini(n,eq%ceq)
  end subroutine pytqini
  
  
  subroutine pytqrfil(filename,eq)
    implicit none
    character(*) :: filename
    type(eq_wrapper) :: eq
    call tqrfil(filename,eq%ceq)
  end subroutine pytqrfil
  
  subroutine pytqrpfil(filename,nsel,selel,eq)
    implicit none
    character(*) :: filename
    integer :: nsel
    character(2) :: selel(:)
    type(eq_wrapper) :: eq
    call tqrpfil(filename,nsel,selel,eq%ceq)
  end subroutine pytqrpfil
  
  subroutine pytqgcom(comp,eq)
    implicit none
    type(comp_wrapper), intent(out) :: comp
    type(eq_wrapper) :: eq
    call tqgcom(comp%n,comp%compnames,eq%ceq)
  end subroutine pytqgcom
  
  subroutine pytqgnp(n,eq)
    implicit none
    integer, intent(out) :: n
    type(eq_wrapper) :: eq
    call tqgnp(n,eq%ceq)
  end subroutine pytqgnp
  
  subroutine pytqgpn(phtupx,phasename,eq)
    implicit none
    integer :: phtupx
    character(*), intent(out) :: phasename
    type(eq_wrapper) :: eq
    call tqgpn(phtupx,phasename,eq%ceq)
  end subroutine pytqgpn
    
  subroutine pytqgpi(phtupx,phasename,eq)
    implicit none
    integer, intent(out) :: phtupx
    character(*) :: phasename
    type(eq_wrapper) :: eq
    call tqgpi(phtupx,phasename,eq%ceq)
  end subroutine pytqgpi
  
  subroutine pytqgpi2(iph,ics,phasename,eq)
    implicit none
    integer, intent(out)  :: iph, ics
    character(*) :: phasename
    type(eq_wrapper) :: eq
    call tqgpi2(iph,ics,phasename,eq%ceq)
  end subroutine pytqgpi2
    
  subroutine pytqgpcn2(n,c,csname)
    implicit none
    integer :: n
    integer :: c
    character(*), intent(out) :: csname
    call tqgpcn2(n,c,csname)
  end subroutine pytqgpcn2
  
  subroutine pytqgpcs(c,nspel,ielno,stoi,smass,qsp)
    implicit none
    integer :: c
    integer, intent(out) :: nspel
    integer :: ielno(*)
    double precision :: stoi(*)
    double precision, intent(out) :: smass, qsp
    call tqgpcs(c,nspel,ielno,stoi,smass,qsp)
  end subroutine pytqgpcs
  
  subroutine pytqphsts(phtupx,newstat,val,eq)
    integer :: phtupx,newstat
    double precision :: val
    type(eq_wrapper) :: eq
    call tqphsts(phtupx,newstat,val,eq%ceq)
  end subroutine pytqphsts
  
  subroutine pytqphsts2(phnames,newstat,val,eq)
    character(*) :: phnames
    integer :: newstat
    double precision :: val
    type(eq_wrapper) :: eq
    call tqphsts2(phnames,newstat,val,eq%ceq)
  end subroutine pytqphsts2
    
  subroutine pytqsetc(stavar,nn1,nn2,value,cnum,eq)
    implicit none
    integer :: nn1
    integer :: nn2
    integer, intent(out) :: cnum
    character(*) :: stavar
    double precision :: value
    type(eq_wrapper) :: eq
    call tqsetc(stavar,nn1,nn2,value,cnum,eq%ceq)
  end subroutine pytqsetc
  
  subroutine pytqtgsw(i)
    implicit none
    integer :: i
    call tqtgsw(i)
  end subroutine pytqtgsw
  
  subroutine pytqce(target,nn1,nn2,value,eq)
    implicit none
    integer :: nn1,nn2
    character(*) :: target
    double precision value
    type(eq_wrapper) :: eq
    call tqce(target,nn1,nn2,value,eq%ceq)
  end subroutine pytqce
  
  subroutine pytqgetv(stavar,nn1,nn2,nn3in,nn3out,values,eq)
    implicit none
    integer :: nn1,nn2,nn3in
    integer, intent(out) :: nn3out
    character(*) ::  stavar
    double precision :: values(*)
    type(eq_wrapper) :: eq
    nn3out=nn3in
    call tqgetv(stavar,nn1,nn2,nn3out,values,eq%ceq)
  end subroutine pytqgetv
  
  subroutine pytqgphc1(iph,nsub,cinsub,spix,yfrac,sites,extra,eq)
    implicit none
    integer :: iph
    integer, intent(out) :: nsub
    integer :: cinsub(*),spix(*)
    double precision :: sites(*),yfrac(*),extra(*)
    type(eq_wrapper) :: eq
    call tqgphc1(iph,nsub,cinsub,spix,yfrac,sites,extra,eq%ceq)
  end subroutine pytqgphc1
  
  subroutine pytqsphc1(nn1,yfra,extra,eq)
    implicit none
    integer :: nn1
    double precision :: yfra(*),extra(*)
    type(eq_wrapper) :: eq
    call tqsphc1(nn1,yfra,extra,eq%ceq)
  end subroutine pytqsphc1
  
  subroutine pytqcph1(nn1,nn2,nn3,gtp,dgdy,d2gdydt,d2gdydp,d2gdy2,eq)
    implicit none
    integer :: nn1,nn2,nn3
    double precision :: gtp(6),dgdy(*),d2gdydt(*),d2gdydp(*),d2gdy2(*)
    type(eq_wrapper) :: eq
    call tqcph1(nn1,nn2,nn3,gtp,dgdy,d2gdydt,d2gdydp,d2gdy2,eq%ceq)
  end subroutine pytqcph1
  
  subroutine pytqcph2(nn1,nn2,nn3,nn4,eq)
    implicit none
    integer :: nn1,nn2,nn3,nn4
    type(eq_wrapper) :: eq
    call tqcph2(nn1,nn2,nn3,nn4,eq%ceq)
  endsubroutine pytqcph2
  
  subroutine pytqcph3(nn1,nn2,g,eq)
    implicit none
    integer :: nn1,nn2
    double precision :: g(*)
    type(eq_wrapper) :: eq
    call tqcph3(nn1,nn2,g,eq%ceq)
  endsubroutine pytqcph3
  
  subroutine pytqdceq(name)
    implicit none
    character(24) :: name
    call tqdceq(name)
  endsubroutine pytqdceq
  
  subroutine pytqcceq(name,nn1,neweq,eq)
    implicit none
    character(24) :: name
    integer, intent(out) :: nn1
    type(eq_wrapper), intent(out) :: neweq
    type(eq_wrapper) :: eq
    call tqcceq(name,nn1,neweq%ceq,eq%ceq)
  endsubroutine pytqcceq
  
  subroutine pytqselceq(name,eq)
    implicit none
    character(24) :: name
    type(eq_wrapper), intent(out) :: eq
    call tqselceq(name,eq%ceq)
  end subroutine pytqselceq
  
  subroutine pytqcref(ciel,phase,tpref,eq)
    implicit none
    integer :: ciel
    character(*) :: phase
    double precision :: tpref(*)
    type(eq_wrapper) :: eq
    call tqcref(ciel,phase,tpref,eq%ceq)
  end subroutine pytqcref
  
  subroutine pytqlr(lut,eq)
    implicit none
    integer :: lut
    type(eq_wrapper) :: eq
    call tqlr(lut,eq%ceq)
  end subroutine pytqlr
  
  subroutine pytqlc(lut,eq)
    implicit none
    integer :: lut
    type(eq_wrapper) :: eq
    call tqlc(lut,eq%ceq)
  end subroutine pytqlc
  
  subroutine pytqquiet(yes)
    implicit none
    logical :: yes
    call tqquiet(yes)
  end subroutine pytqquiet
  
end module RawOpenCalphad
