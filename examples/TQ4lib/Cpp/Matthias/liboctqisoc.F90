module cstr

    contains

    function c_to_f_string(s) result(str)
        use iso_c_binding
        implicit none
        character(kind=c_char,len=1), intent(in) :: s(*)
        character(len=:), allocatable :: str
        integer i, nchars
        i = 1
        do
            if (s(i) == c_null_char) exit
            i = i + 1
        end do
        nchars = i - 1
        allocate(character(len=nchars) :: str)
        str = transfer(s(1:nchars), str)
    end function c_to_f_string

    subroutine f_to_c_string(fstring, cstr)
        use iso_c_binding
        implicit none
        character(len=24) :: fstring
        character(kind=c_char, len=1), intent(out) :: cstr(*)
        integer i
        do i = 1, len(fstring)
            cstr(i) = fstring(i:i)
            cstr(i+1) = c_null_char
        end do
    end subroutine f_to_c_string

end module cstr

module liboctqisoc
    use iso_c_binding
    use cstr
    use liboctq
    implicit none
    integer(c_int), bind(c) :: c_nel
    integer(c_int), bind(c) :: c_maxc=20
    integer(c_int), bind(c) :: c_maxp=100
    type(c_ptr), bind(c), dimension(maxc) :: c_cnam
    character(len=25), dimension(maxc), target :: cnames
    integer(c_int), bind(c) :: c_ntup

    TYPE, bind(c) :: c_gtp_equilibrium_data
        integer(c_int) :: status,multiuse,eqno,next
        character(c_char) :: eqname*24
        character(c_char) :: comment*72
        real(c_double) :: tpval(2)
        real(c_double) :: rtn
        real(c_double) :: weight
        real(c_double) :: svfunres
        TYPE(c_ptr) :: lastcondition
        TYPE(c_ptr) :: lastexperiment
        TYPE(c_ptr) :: complist
        real(c_double) :: compstoi
        real(c_double) :: invcompstoi
        TYPE(c_ptr) :: phase_varres
        TYPE(c_ptr) :: eq_tpres
        real(c_double) :: cmuval
        real(c_double) :: xconv
        real(c_double) :: gmindif=-5.0D-2
        integer(c_int) :: maxiter
        character(c_char) :: eqextra*80
        integer(c_int) :: sysmatdim=0
        integer(c_int) :: nfixmu=0
        integer(c_int) :: nfixph=0
        integer(c_int) :: fixmu
        integer(c_int) :: fixph
        real(c_double) :: savesysmat
    END TYPE c_gtp_equilibrium_data

    contains

    integer function c_noofcs(iph) bind(c, name='c_noofcs')
        integer(c_int), value :: iph
        c_noofcs = noofcs(iph)
        return
    end function c_noofcs

    integer function c_ierr() bind(c, name='c_ierr')
        c_ierr=gx%bmperr
        return
    end function c_ierr

    subroutine c_tqini(n, c_ceq) bind(c, name='c_tqini')
        integer(c_int), intent(in) :: n
        type(c_ptr), intent(out) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        !=================
        call tqini(n, ceq)
        !=================
        c_ceq = c_loc(ceq)
    end subroutine c_tqini

    subroutine c_tqrfil(filename,c_ceq) bind(c, name='c_tqrfil')
        character(kind=c_char,len=1), intent(in) :: filename(*)
        character(len=:), allocatable :: fstring
        type(gtp_equilibrium_data), pointer :: ceq
        type(c_ptr), intent(inout) :: c_ceq
        integer :: i
        integer :: j
        integer :: l
        character(kind=c_char, len=1),dimension(24), target :: f_pointers
        call c_f_pointer(c_ceq, ceq)
        fstring = c_to_f_string(filename)
        !========================
        call tqrfil(fstring, ceq)
        !========================
        c_ntup = ntup
        c_nel = nel
        do i = 1, nel
            cnames(i) = trim(cnam(i)) // c_null_char
            c_cnam(i) = c_loc(cnames(i))
        end do
        c_ceq = c_loc(ceq)
    end subroutine c_tqrfil

    subroutine c_tqrpfil(filename,nel,c_selel,c_ceq) bind(c, name='c_tqrpfil')
        character(kind=c_char), intent(in) :: filename
        integer(c_int), intent(in), value :: nel
        type(c_ptr), intent(in), dimension(nel), target :: c_selel
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        character(len=:), allocatable :: fstring
        character, pointer :: selel(:)
        integer :: i
        character elem(nel)*2
        fstring = c_to_f_string(filename)
        call c_f_pointer(c_ceq, ceq)
        do i = 1, nel
            call c_f_pointer(c_selel(i), selel, [3])
            elem(i) = c_to_f_string(selel)
        end do
        !====================================
        call tqrpfil(fstring, nel, elem, ceq)
        !====================================
        c_ntup = ntup
        c_nel = nel
        do i = 1, nel
            cnames(i) = trim(cnam(i)) // c_null_char
            c_cnam(i) = c_loc(cnames(i))
        end do
        c_ceq = c_loc(ceq)
    end subroutine c_tqrpfil

    subroutine c_tqgcom(n,components,c_ceq) bind(c, name='c_tqgcom')
        integer(c_int), intent(inout) :: n
        type(c_ptr), intent(inout) :: c_ceq
        character(kind=c_char, len=1), intent(out) :: components(maxel*3)
        integer, target :: nc
        character(len=24) :: fcomponents(maxel)
        type(gtp_equilibrium_data), pointer :: ceq
        integer :: i,j,l
        call c_f_pointer(c_ceq, ceq)
        !================================
        call tqgcom(nc, fcomponents, ceq)
        !================================
        l = 1
        do i = 1, nc
            do j = 1, 2
                components(l)(1:1) = fcomponents(i)(j:j)
                l=l+1
            end do
        end do
        ! null termination
        components(i*2-1) = c_null_char
        c_ceq = c_loc(ceq)
        n = nc
    end subroutine c_tqgcom

    subroutine c_tqgnp(n, c_ceq) bind(c, name='c_tqgnp')
        integer(c_int), intent(inout) :: n
        type(c_ptr), intent(inout) :: c_ceq
        integer, target :: nc
        type(gtp_equilibrium_data), pointer :: ceq
        call c_f_pointer(c_ceq, ceq)
        !=================
        call tqgnp(n, ceq)
        !=================
        c_ceq = c_loc(ceq)

    end subroutine c_tqgnp

    subroutine c_tqgpn(n,phasename, c_ceq) bind(c, name='c_tqgpn')
        integer(c_int), intent(in), value :: n
        character(kind=c_char, len=1), intent(inout) :: phasename(36)
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        character(len=24) :: fstring
        integer :: i
        call c_f_pointer(c_ceq, ceq)
        !==========================
        call tqgpn(n, fstring, ceq)
        !==========================
        do i=1,len(trim(fstring))
            phasename(i)(1:1) = fstring(i:i)
            phasename(i+1)(1:1) = c_null_char
        end do
        c_ceq = c_loc(ceq)
    end subroutine c_tqgpn

    subroutine c_tqgpi(n,phasename,c_ceq) bind(c, name='c_tqgpi')
        integer(c_int), intent(out) :: n
        character(c_char), intent(in) :: phasename(24)
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        character(len=24) :: fstring
        call c_f_pointer(c_ceq, ceq)
        fstring = c_to_f_string(phasename)
        !==========================
        call tqgpi(n, fstring, ceq)
        !==========================
        c_ceq = c_loc(ceq)
    end subroutine c_tqgpi

    subroutine c_tqgpcn2(n, c, constituentname, c_ceq) bind(c, name='c_tqgpcn2')
        integer(c_int), intent(in), value :: n
        integer(c_int), intent(in), value :: c
        character(c_char), intent(out) :: constituentname(24)
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        character fstring*(24)
        double precision mass
        call c_f_pointer(c_ceq, ceq)
        !==================================================
        call get_constituent_name(n,c,fstring,mass)
        !==================================================
        call f_to_c_string(fstring, constituentname)
        c_ceq = c_loc(ceq)
    end subroutine c_tqgpcn2

    subroutine c_tqgpcn(n, c, constituentname, c_ceq) bind(c, name='c_tqgpcn')
        integer(c_int), intent(in), value :: n
        integer(c_int), intent(in), value :: c
        character(c_char), intent(out) :: constituentname(24)
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        character fstring*(24)
        double precision mass
        integer :: i
        call c_f_pointer(c_ceq, ceq)
        !==========================================
        call get_constituent_name(n,c,fstring,mass)
        !==========================================
        !call f_to_c_string(fstring, constituentname)
        do i=1,len(trim(fstring))
            constituentname(i)(1:1) = fstring(i:i)
            constituentname(i+1)(1:1) = c_null_char
        end do
        c_ceq = c_loc(ceq)
    end subroutine c_tqgpcn

    subroutine c_tqgpci(n,c, constituentname, c_ceq) bind(c, name='c_tqgpci')
        integer(c_int), intent(in) :: n
        integer(c_int), intent(out) :: c
        character(c_char), intent(in) :: constituentname(24)
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        character(len=24) :: fstring
        fstring = c_to_f_string(constituentname)
        call c_f_pointer(c_ceq, ceq)
        !==============================
        call tqgpci(n, c, fstring, ceq)
        !==============================
        c_ceq = c_loc(ceq)
    end subroutine c_tqgpci

    subroutine c_tqgpcs(n, c, stoi, mass, c_ceq) bind(c, name='c_tqgpcs')
        integer(c_int), intent(in) :: n
        integer(c_int), intent(in) :: c
        real(c_double), intent(out) :: stoi(*)
        real(c_double), intent(out) :: mass
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        call c_f_pointer(c_ceq, ceq)
        !=============================
        call tqgpcs(n,c,stoi,mass,ceq)
        !=============================
        c_ceq=c_loc(ceq)
    end subroutine c_tqgpcs

    subroutine c_tqgccf(n1,n2,elnames,stoi,mass,c_ceq) bind(c, name='c_tqgccf')
        integer(c_int), intent(in) :: n1
        integer(c_int), intent(out) :: n2
        character(c_char), intent(out) :: elnames(2)
        real(c_double), intent(out) :: stoi(*)
        real(c_double), intent(out) :: mass
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        call c_f_pointer(c_ceq, ceq)
        !=========================================
        call tqgccf(n1,n2,elnames,stoi, mass, ceq)
        !=========================================
        c_ceq = c_loc(ceq)
    end subroutine c_tqgccf

    subroutine c_tqgnpc(n,c,c_ceq) bind(c, name='c_tqgnpc')
        integer(c_int), intent(in) :: n
        integer(c_int), intent(out) :: c
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        call c_f_pointer(c_ceq,ceq)
        !===================
        call tqgnpc(n,c,ceq)
        !===================
        c_ceq = c_loc(ceq)
    end subroutine c_tqgnpc

    subroutine c_tqphtupsts(phtupx,newstat,val,c_ceq) &
            bind(c, name='c_tqphtupsts')
        integer(c_int), intent(in), value :: phtupx
        integer(c_int), intent(in), value :: newstat
        real(c_double), intent(in), value :: val
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        call c_f_pointer(c_ceq,ceq)
        !======================================
        call tqphtupsts(phtupx,newstat,val,ceq)
        !======================================
        c_ceq = c_loc(ceq)
    end subroutine c_tqphtupsts

    subroutine c_tqsetc(statvar, n1, n2, mvalue, cnum, c_ceq) &
            bind(c, name='c_tqsetc')
        integer(c_int), intent(in),value :: n1
        integer(c_int), intent(in),value :: n2 !
        integer(c_int), intent(out) :: cnum
        character(c_char), intent(in) :: statvar
        real(c_double), intent(in), value :: mvalue
        type(gtp_equilibrium_data), pointer :: ceq
        type(c_ptr), intent(inout) :: c_ceq
        call c_f_pointer(c_ceq, ceq)
        !==============================================
        call tqsetc(statvar, n1, n2, mvalue, cnum, ceq)
        !==============================================
        c_ceq = c_loc(ceq)
    end subroutine c_tqsetc

    subroutine c_tqce(mtarget,n1,n2,mvalue,c_ceq) bind(c,name='c_tqce')
        integer(c_int), intent(in),value :: n1
        integer(c_int), intent(in),value :: n2
        type(c_ptr), intent(inout) :: c_ceq
        character(c_char), intent(inout) :: mtarget
        real(c_double), intent(inout) :: mvalue
        type(gtp_equilibrium_data), pointer :: ceq
        character(len=24) :: fstring
        call c_f_pointer(c_ceq,ceq)
        fstring = c_to_f_string(mtarget)
        !==================================
        call tqce(fstring,n1,n2,mvalue,ceq)
        !==================================
        c_ceq = c_loc(ceq)
    end subroutine c_tqce

    subroutine c_tqgetv(statvar,n1,n2,n3,values,c_ceq) bind(c,name='c_tqgetv')
        integer(c_int), intent(in), value ::  n1
        integer(c_int), intent(in), value ::  n2
        integer(c_int), intent(inout) :: n3
        character(c_char), intent(in) :: statvar
        real(c_double), intent(inout) :: values(*)
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        character(len=24) :: fstring
        !integer :: n
        integer :: i
        call c_f_pointer(c_ceq, ceq)
        fstring = c_to_f_string(statvar)
        !========================================
        call tqgetv(fstring, n1,n2,n3,values,ceq)
        !========================================
        c_ceq = c_loc(ceq)
    end subroutine c_tqgetv

    subroutine c_tqgphc1(n1,nsub,cinsub,spix,yfrac,sites,extra,c_ceq) &
            bind(c,name='c_tqgphc1')
        integer(c_int), intent(in), value :: n1
        integer(c_int), intent(out) :: nsub
        integer(c_int), intent(out) :: cinsub(*)
        integer(c_int), intent(in) :: spix(*)
        real(c_double), intent(in) :: sites(*)
        real(c_double), intent(in) :: yfrac(*)
        real(c_double), intent(in) :: extra(*)
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        call c_f_pointer(c_ceq, ceq)
        !======================================================
        call tqgphc1(n1,nsub,cinsub,spix,yfrac,sites,extra,ceq)
        !======================================================
        c_ceq = c_loc(ceq)
    end subroutine c_tqgphc1

    subroutine c_tqsphc1(n1,yfra,extra,c_ceq) bind(c,name='c_tqsphc1')
        integer(c_int), intent(in), value :: n1
        real(c_double), intent(in) ::yfra(*)
        real(c_double), intent(out) :: extra(*)
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        call c_f_pointer(c_ceq, ceq)
        !=====================================================================
        call set_constitution(phasetuple(n1)%ixphase,phasetuple(n1)%compset, &
                              yfra,extra,ceq)
        !=====================================================================
        c_ceq = c_loc(ceq)
    end subroutine c_tqsphc1

    subroutine c_tqcph1(n1,n2,n3,gtp,dgdy,d2gdydt,d2gdydp,d2gdy2,c_ceq) &
            bind(c,name='c_tqcph1')
        integer(c_int), intent(in), value :: n1
        integer(c_int), intent(in), value :: n2
        integer(c_int), intent(out) :: n3
        real(c_double), intent(out) :: gtp(6)
        real(c_double), intent(out) :: dgdy(*)
        real(c_double), intent(out) :: d2gdydt(*)
        real(c_double), intent(out) :: d2gdydp(*)
        real(c_double), intent(out) :: d2gdy2(*)
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        call c_f_pointer(c_ceq, ceq)
        !========================================================
        call tqcph1(n1,n2,n3,gtp,dgdy,d2gdydt,d2gdydp,d2gdy2,ceq)
        !========================================================
        c_ceq = c_loc(ceq)
    end subroutine c_tqcph1

    subroutine c_tqcph2(n1,n2,n3,n4,c_ceq) bind(c,name='c_tqcph2')
        integer(c_int), intent(in), value :: n1
        integer(c_int), intent(in), value :: n2
        integer(c_int), intent(out) :: n3
        integer(c_int), intent(out) :: n4
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        call c_f_pointer(c_ceq, ceq)
        !===========================
        call tqcph2(n1,n2,n3,n4,ceq)
        !===========================
        c_ceq = c_loc(ceq)
    end subroutine c_tqcph2

    subroutine c_tqdceq(ceqname) bind(c,name='c_tqdceq')
        character(c_char), intent(in) :: ceqname(24)
        character(len=24) :: fstring
        fstring = c_to_f_string(ceqname)
        !================
        call tqdceq(fstring)
        !================
    end subroutine c_tqdceq

    subroutine c_tqcceq(ceqname,n1,c_newceq,c_ceq) bind(c,name='c_tqcceq')
        character(c_char), intent(in) :: ceqname(24)
        integer(c_int), intent(out) :: n1
        type(c_ptr), intent(inout) :: c_newceq
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        type(gtp_equilibrium_data), pointer :: newceq
        character(len=24) :: fstring
        call c_f_pointer(c_newceq, newceq)
        call c_f_pointer(c_ceq, ceq)
        fstring = c_to_f_string(ceqname)
        !=================================
        call tqcceq(fstring,n1,newceq,ceq)
        !=================================
        c_newceq = c_loc(newceq)
        c_ceq = c_loc(ceq)
    end subroutine c_tqcceq

    subroutine c_tqselceq(ceqname,c_ceq) bind(c,name='c_tqselceq')
        character(c_char), intent(in) :: ceqname(24)
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        character(len=24) :: fstring
        call c_f_pointer(c_ceq, ceq)
        fstring = c_to_f_string(ceqname)
        !=========================
        call tqselceq(fstring,ceq)
        !=========================
        c_ceq = c_loc(ceq)
    end subroutine c_tqselceq

    subroutine c_reset_conditions(cline,c_ceq) bind(c,name='c_reset_conditions')
        character(c_char), intent(in) :: cline(24)
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        character(len=24) :: fstring
        call c_f_pointer(c_ceq, ceq)
        fstring = c_to_f_string(cline)
        !=================================
        call reset_conditions(fstring,ceq)
        !=================================
        c_ceq = c_loc(ceq)
    end subroutine c_reset_conditions

    subroutine c_Change_Status_Phase(myname,nystat,myval,c_ceq) &
            bind(c,name='c_Change_Status_Phase')
        character(c_char), intent(in) :: myname(24)
        integer(c_int), intent(in), value :: nystat
        real(c_double), intent(in), value :: myval
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        character(len=24) :: fstring
        call c_f_pointer(c_ceq, ceq)
        fstring = c_to_f_string(myname)
        !=================================================
        call Change_Status_Phase(fstring,nystat,myval,ceq)
        !=================================================
        c_ceq = c_loc(ceq)
    end subroutine c_Change_Status_Phase

    subroutine c_tqlr(lut,c_ceq) bind(c,name='c_tqlr')
        integer(c_int), intent(in), value :: lut
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        call c_f_pointer(c_ceq, ceq)
        !=================
        call tqlr(lut,ceq)
        !=================
        c_ceq = c_loc(ceq)
    end subroutine c_tqlr

    subroutine c_tqlc(lut,c_ceq) bind(c,name='c_tqlc')
        integer(c_int), intent(in), value :: lut
        type(c_ptr), intent(inout) :: c_ceq
        type(gtp_equilibrium_data), pointer :: ceq
        call c_f_pointer(c_ceq, ceq)
        !=================
        call tqlc(lut,ceq)
        !=================
        c_ceq = c_loc(ceq)
    end subroutine c_tqlc

end module liboctqisoc
