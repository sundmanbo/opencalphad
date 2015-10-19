!
! gtp3Z included in gtp3.F90
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!***************************************************************
! library with TP functions used by general thermodynamic package
!
! the declarations below are all moved to gtp3Z.F90
!
! MODULE TPFUNLIB
!
! Copyright 2009-2015, Bo Sundman, France
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!-------------------------------------------------------------------------
!
!  use METLIB
!
!  implicit double precision (a-h,o-z)
! error codes, copied from GTP, currently not used
!  character (len=64), dimension(4000:4029) :: tpferrmess=&
!      ['Too many coefficients in a TP function.                         ',&
!       'Illegal character in a TP function, digit expected.             ',&
!       'Unkown symbol in TP function                                    ',&
!       'Expected ( after unary function                                 ',&
!       'Too many ) in a TP function                                     ',&
!       'Illegal character in a TP function                              ',&
!       'Too few ) in a TP function                                      ',&
!       'Too many ( in exponent                                          ',&
!       'Illegally placed ( in the exponent of a TP function             ',&
!       'No digits after ( in the exponent of a TP function              ',&
! 4010:
!       'Illegally placed ) in exponent in TP function                   ',&
!       'Too high power in a TP function, max 99, min -99                ',&
!       'Missing ) in the exponent of a TP function                      ',&
!       'Illegal termination of a TP function reading a TDB file         ',&
!       'No more free TP root records                                    ',&
!       'No more free TP expression records                              ',&
!       'Illegal expression inside unary argument of a TP function       ',&
!       'Illegal code found when evaluating a TP function                ',&
!       'Found acoefficent zero in a term of a TP function               ',&
!       'Illegal code in a TP function                                   ',&
! 4020:
!       'Negative argument to logarithm in a TP function                 ',&
!       'Unknown unary function in evaluation for a TP function          ',&
!       'Too many symbols in a TP function term                          ',&
!       'Two unary functions in a TP function term                       ',&
!       'Too complicated TP function term                                ',&
!       'Too many temperature ranges in a TP function                    ',&
!       'TP function with same name already entered                      ',&
!       'Symbol referenced in a parameter does not exist                 ',&
!       'Missing separator between phase and constituent array in paramet',&
!       'Cannot enter disordered fraction set when several composition se']
!
! length of a function symbol
!  PARAMETER (lenfnsym=16)
!  integer, private :: freetpfun
!  PARAMETER (zero=0.0D0,one=1.0D0)
!-----------------------------------------------------------------
!-\begin{verbatim}
!  TYPE gtp_parerr
! This record contains the global error code.  In parallell processing each
! parallell processes has its own error code copied to this if nonzero
! it should be replaced by gtperr for separate errors in treads
!     INTEGER :: bmperr
!  END TYPE gtp_parerr
!  TYPE(gtp_parerr) :: gx
! needed to have error code as private in threads
!-$OMP  threadprivate(gx)
!-\end{verbatim}
!-----------------------------------------------------------------
!-\begin{verbatim}
!  integer, parameter :: tpfun_expression_version=1
!  TYPE tpfun_expression
! Coefficients, T and P powers, unary functions and links to other functions
!     integer noofcoeffs,nextfrex
!     double precision, dimension(:), pointer :: coeffs
! each coefficient kan have powers of T and P/V and links to other TPFUNS
! and be multiplied with a following LOG or EXP term. 
!     integer, dimension(:), pointer :: tpow
!     integer, dimension(:), pointer :: ppow
!     integer, dimension(:), pointer :: wpow
!     integer, dimension(:), pointer :: plevel
!     integer, dimension(:), pointer :: link
!  END TYPE tpfun_expression
! These records are allocated when needed, not stored in arrays
!-\end{verbatim}
!-----------------------------------------------------------------
! BITS in TPFUN
! TPCONST     set if a constant value
! TPOPTCON    set if optimizing value
! TPNOTENT    set if referenced but not entered (when reading TDB files)
! TPVALUE     set if evaluated only explicitly (keeping its value)
!  integer, parameter :: &
!       TPCONST=0,    TPOPTCON=1,   TPNOTENT=2,    TPVALUE=3
!-----------------------------------------------------------------
!-\begin{verbatim}
!  integer, parameter :: tpfun_root_version=1
!  TYPE tpfun_root
! Root of a TP function including name with links to coefficients and codes
! and results.  Note that during calculations which can be parallellized
! the results can be different for each parallell process
!     character*(lenfnsym) symbol
! limits is the low temperature limit for each range
! funlinks links to expression records for each range
! each range can have its own function, status indicate if T and P or T and V
!     integer noofranges,nextfree,status
!     double precision, dimension(:), pointer :: limits
!     TYPE(tpfun_expression), dimension(:), pointer :: funlinks
!     double precision hightlimit
!  END TYPE tpfun_root
! These records are stored in arrays as the actual function is global but each
! equilibrium has its own result array (tpfun_parres) depending on the local
! values of T and P/V.  The same indiex is used in the global and local arrays.
! allocated in init_gtp
!  TYPE(tpfun_root), private, dimension(:), pointer :: tpfuns
!-\end{verbatim}
!
!-----------------------------------------------------------------
!-\begin{verbatim}
!  integer, parameter :: tpfun_parres_version=1
!  TYPE tpfun_parres
! Contains a TP results, 6 double for results and 2 doubles for T and P 
! values used to calculate the results
! Note that during calculations which can be parallellized the final
! results can be different for each parallell process
!     double precision, dimension(2) :: tpused
!     double precision, dimension(6) :: results
!  END TYPE tpfun_parres
! This array is local to the gtp_equilibrium_data record
! index the same as the function
!-\end{verbatim}
!
!-----------------------------------------------------------------
!
!CONTAINS
!
!\begin{verbatim}
 SUBROUTINE tpfun_init(nf,tpres)
! allocate tpfuns and create a free list inside the tpfuns
   implicit none
   integer nf
! use tpres declared externally for parallel processing
  TYPE(tpfun_parres), dimension(:), pointer :: tpres
!\end{verbatim}
   integer ifri
   allocate(tpfuns(nf))
   allocate(tpres(nf))
! create free list for named functions records
   freetpfun=1
   do ifri=1,nf-1
      tpfuns(ifri)%nextfree=ifri+1
      tpfuns(ifri)%noofranges=0
   enddo
   tpfuns(nf)%nextfree=-1
   return
 END SUBROUTINE tpfun_init

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 integer function notpf()
! number of tpfunctions because freetpfun is private
   implicit none
!\end{verbatim}
   notpf=freetpfun-1
 end function notpf

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  SUBROUTINE find_tpfun_by_name(name,lrot)
! returns the location of a TP function
! if lrot>0 then start after lrot, this is to allow finding with wildcard *
    implicit none
    integer lrot
    character name*(*)
!\end{verbatim} %+
    character name1*16
    integer i,j
    name1=name
    call capson(name1)
    if(lrot.le.0 .or. lrot.ge.freetpfun) then
       j=1
    else
! if 1 < lrot < freetpfun start looking from lrot+1
       j=lrot+1
    endif
    do i=j,freetpfun-1
       if(compare_abbrev(name,tpfuns(i)%symbol)) then
          lrot=i; goto 1000
       endif
    enddo
    gx%bmperr=4060
1000 continue
    return
  end SUBROUTINE find_tpfun_by_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
  SUBROUTINE find_tpfun_by_name_exact(name,lrot,notent)
! returns the location of a TP function, notent TRUE if not entered
    implicit none
    integer lrot
    logical notent
    character name*(*)
!\end{verbatim}
    character name1*16
    integer i
    notent=.FALSE.
    name1=name
    call capson(name1)
    do i=1,freetpfun-1
       if(name.eq.tpfuns(i)%symbol) then
          lrot=i
          if(btest(tpfuns(i)%status,TPNOTENT)) then
             notent=.TRUE.
          endif
          goto 1000
       endif
    enddo
    gx%bmperr=4060
1000 continue
    return
  end SUBROUTINE find_tpfun_by_name_exact

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine eval_tpfun(lrot,tpval,result,tpres)
!    subroutine eval_tpfun(lrot,tpval,symval,result)
! evaluate a TP function with several T ranges
!   implicit double precision (a-h,o-z)
   implicit none
   integer lrot
   double precision tpval(2),result(6)
   TYPE(tpfun_parres), dimension(:), pointer :: tpres
!\end{verbatim}
   integer nr,ns
   TYPE(tpfun_expression), pointer :: exprot
! mini is the maximum relative difference between calculated and current values
! of T and P for using the stored values of a function
   double precision, parameter :: mini=1.0D-8
! use lowest range for all T values lower than first upper limit
! and highest range for all T values higher than the next highest limit
! one should signal if T is lower than lowest limit or higher than highest
! used  saved reults if same T and P
!
! SURPRICE: removing the do loops reduces CPU time with 0.1 seconds ...
   if(lrot.le.0) then
      result=zero
      goto 1000
   elseif(btest(tpfuns(lrot)%status,TPCONST)) then
! TP symbol is a constant
      result=zero
      result(1)=tpfuns(lrot)%limits(1)
      goto 1000
   else
! check if previous values can be used
      if(abs(tpres(lrot)%tpused(1)-tpval(1)).le.mini*tpres(lrot)%tpused(1) &
 .and.(abs(tpres(lrot)%tpused(2)-tpval(2)).le.mini*tpres(lrot)%tpused(2))) &
           then
         result=tpres(lrot)%results
         goto 1000
      else
         result=zero
      endif
   endif
! we must calculate the function
   nr=tpfuns(lrot)%noofranges
   if(nr.eq.1) then
      exprot=>tpfuns(lrot)%funlinks(1)
      call ct1efn(exprot,tpval,result,tpres)
   else
      ns=1
      do while(ns.lt.nr)
         if(tpval(1).lt.tpfuns(lrot)%limits(ns+1)) then
            exprot=>tpfuns(lrot)%funlinks(ns)
            call ct1efn(exprot,tpval,result,tpres)
! for debug output below
            nr=ns
            goto 900
         endif
         ns=ns+1
      enddo
      exprot=>tpfuns(lrot)%funlinks(nr)
      call ct1efn(exprot,tpval,result,tpres)
   endif
! save the calculated results
900 continue
   tpres(lrot)%tpused(1)=tpval(1)
   tpres(lrot)%tpused(2)=tpval(2)
!   new: do i=1,6
!      tpres(lrot)%results(i)=result(i)
!   enddo new
   tpres(lrot)%results=result
22  format(A,i3,4(1PE11.2))
1000 continue
   return
 end subroutine eval_tpfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_tpfun(lrot,nosym,str)
! lists a TP symbols with several ranges into string str
! lrot is index of function, if nosym=0 the function name is copied to str
   implicit none
   character str*(*)
   integer nosym,lrot
!\end{verbatim} %+
   integer ip,nr
   character line*2048,tps(2)*1
   TYPE(tpfun_expression), pointer :: exprot
! Handle variables
   if(lrot.le.0) then
      str=' =0; N '
      goto 1000
   elseif(btest(tpfuns(lrot)%status,TPCONST)) then
! UNFINISHED temporarily list all optimizing variables
      if(btest(tpfuns(lrot)%status,TPOPTCON)) then
         if(tpfuns(lrot)%limits(1).eq.zero) then
! this is a clumsy way to suppress listing optimizing coeff that are zero
            str='_A00 '; goto 1000
         endif
      endif
      line=tpfuns(lrot)%symbol
      ip=len_trim(line)
      line(ip+1:ip+3)=' = '
      ip=ip+4
      call wrinum(line,ip,12,0,tpfuns(lrot)%limits(1))
      goto 900
   endif
! these are the symbols used to represent T and P
   tps(1)='T'
   tps(2)='P'
   if(nosym.eq.0) then
      line=tpfuns(lrot)%symbol
      ip=len_trim(line)
      line(ip+1:ip+3)=' = '
      ip=ip+4
   else
      line='= '
      ip=3
   endif
   if(lrot.le.0) then
      line(ip:)=' 298.15  0; 6000 N'
      goto 900
   endif
!   nr=1
   do nr=1,tpfuns(lrot)%noofranges
!      write(line(ip:ip+10),10)tpfuns(lrot)%limits(nr)
!10     format(F8.2,' Y ')
!      ip=ip+9
!      write(*,*)'tpfun4: ',lrot,tpfuns(lrot)%noofranges,nr
      if(tpfuns(lrot)%limits(nr).gt.1.0D3) then
         call wrinum(line,ip,8,0,tpfuns(lrot)%limits(nr))
      else
! problem as 298.15 is written as 298.14999 and
! problem that 1 is written as 1.00001
!         call wrinum(line,ip,6,0,tpfuns(lrot)%limits(nr)+1.0D-5)
         call wrinum(line,ip,6,0,tpfuns(lrot)%limits(nr))
      endif
      line(ip:ip+2)=' Y '
      ip=ip+1
      if(nr.gt.1) ip=ip+2
      exprot=>tpfuns(lrot)%funlinks(nr)
      call ct1wfn(exprot,tps,line,ip)
      line(ip:ip+1)='; '
      ip=ip+2
   enddo
!   write(line(ip:ip+10),11)tpfuns(lrot)%hightlimit
!11  format(F8.2,' N ')
!   ip=ip+11
   call wrinum(line,ip,8,0,tpfuns(lrot)%hightlimit)
   line(ip:ip+2)=' N '
900 continue
!   write(*,*)'list_tpfun: ',len(str),len_trim(line)
   if(len_trim(line).gt.len(str)) then
      write(kou,910)' *** WARNING: Character for listing funtion too short',&
           len_trim(line),len(str)
910   format(a,2i5)
   endif
   str=line
20 format(a)
1000 continue
   return
 end subroutine list_tpfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine list_all_funs(lut)
! list all functions except those starting with _ (parameters)
   implicit none
   integer lut
!\end{verbatim}
!   implicit double precision (a-h,o-z)
   integer nosym,ifun
   character str*2048,number*4
   logical once
! nosym=0 means the local symbol name is included in the listing
   once=.TRUE.
   nosym=0
   write(lut,10)
10 format(/'List of all symbols used in phase parameters (TP-functions):'/ &
        ' Predefined symbols:'/&
        ' BELOW(TB) = EXP(20*(1-T/TB))/(1+EXP(20*(1-T/TB)));'/&
        ' ABOVE(TB) = 1-EXP(20*(1-T/TB))/(1+EXP(20*(1-T/TB)));'/&
        ' Nr  Name =     T-low  expression; T-high Y/N')
20  format(I4,1x,A)
!   write(*,*)'First free index: ',freetpfun
   do ifun=1,freetpfun-1
      write(str,20)ifun
      call list_tpfun(ifun,nosym,str(6:))
      if(str(6:9).eq.'_A00 ') then
         if(once) then
            write(lut,30)
30          format(' *** Optimizing coefficents that are zero are not listed')
            once=.FALSE.
         endif
      else
         if(str(6:6).ne.'_') call wrice2(lut,0,12,78,1,str)
      endif
   enddo
   return
 end subroutine list_all_funs

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_unentered_funs(lut,nr)
! counts and list functions with TPNOTENT bit set if lut>0
   implicit none
   integer lut,nr
!\end{verbatim}
!   implicit double precision (a-h,o-z)
   integer nosym,ifun
   nr=0
   do ifun=1,freetpfun-1
      if(btest(tpfuns(ifun)%status,TPNOTENT)) then
         if(lut.gt.0) write(lut,30)tpfuns(ifun)%symbol
30       format('Missing function: ',a)
         nr=nr+1
      endif
   enddo
1000 continue
   return
 end subroutine list_unentered_funs

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 SUBROUTINE ct1xfn(string,ip,nc,coeff,koder,fromtdb)
!...compiles an expression in string from position ip
!     it can refer to T and P or symbols in fnsym
!     compiled expression returned in coeff and koder
!
! >>> this is very messy
!
!...algorithm for function extraction
! 10*T**2 -5*T*LOG(T) +4*EXP(-5*T**(-1))
!
! AT LABEL 100 start of expression or after (
! sign=1
! -, sign=-1                              goto 200
! +, skip
!
! AT LABEL 200 after sign
! if A-Z                                  goto 300
! if 0-9, extract number                  goto 400
! (                                       goto 100
! ;                                       END or ERROR
! empty                                   END or ERROR
! anything else                           ERROR
!
! AT LABEL 300 symbol
! if T or P, extract power if any incl () goto 400
! unary fkn? extract (                    goto 100
! symbol                                  goto 400
!
! AT LABEL 400 after factor
! -, sign=-1                              goto 200
! sign=1
! +, skip                                 goto 200
! )                                       goto 400
! ** or ^ extract and store power incl () goto 400
! *                                       goto 200
! empty                                   goto 900
!
! for TDB compatibility skip #
!
! allow unary functions ABOVE(TB) and BELOW(TB) where TB= is the break temp
! check consistency
   implicit none
   integer ip,nc,koder(5,*)
   character string*(*)
   double precision coeff(*)
   logical fromtdb
!\end{verbatim} %+
!   implicit double precision (a-h,o-z)
   integer, parameter :: nunary=6
   integer, parameter :: lenfnsym=16
   double precision, parameter :: zero=0.0D0,one=1.0D0
   integer i,j,jss,levelp,mterm,ipower,nterm
   double precision sign,val
   character ch1*1
   logical zeroc
   character symbol*(lenfnsym),unary(nunary)*6
   character, parameter :: tsym='T ',psym='P '
   DATA unary/'LOG   ','LN    ','EXP   ','ERF   ','ABOVE ','BELOW '/
!
! coeff(nterm)   double with coefficient
! koder(1,nterm) power of T
! koder(2,nterm) power of P
! koder(3,nterm) power of linked symbol (see koder(5,nterm)
! koder(4,nterm) level of parenthesis
! koder(5,nterm) symbol link or -(unary function index)
   mterm=nc
   levelp=0
   nterm=1
   coeff(1)=zero
   do i=1,5
      koder(i,1)=0
   enddo
!
!...start of expression or after(
100 if(eolch(string,ip)) goto 800
   zeroc=.FALSE.
   ch1=biglet(string(ip:ip))
   sign=one
   if(ch1.eq.'-') then
      sign=-one
      ip=ip+1
      if(coeff(nterm).ne.zero) nterm=nterm+1
      if(nterm.gt.mterm) then
         gx%bmperr=4000
         goto 1000
      endif
      coeff(nterm)=zero
      do i=1,5
         koder(i,nterm)=0
      enddo
   endif
   if(ch1.eq.'+') then
      ip=ip+1
      if(coeff(nterm).ne.zero) nterm=nterm+1
      if(nterm.gt.mterm) then
         gx%bmperr=4000
         goto 1000
      endif
      coeff(nterm)=zero
      do i=1,5
         koder(i,nterm)=0
      enddo
   endif
!
!...allowed: unsigned number or symbol (any previous sign in "sign")
200 continue
   if(eolch(string,ip)) goto 800
   ch1=biglet(string(ip:ip))
   if(ch1.eq.'(') then
      levelp=levelp+1
      if(nterm.eq.0) then
         nterm=1
         coeff(nterm)=zero
         do i=1,5
            koder(i,nterm)=0
         enddo
      endif
      koder(4,nterm)=levelp
      ip=ip+1
      goto 100
   elseif(ch1.eq.';') then
      goto 900
   endif
   if(ch1.ge.'A' .and. ch1.le.'Z') goto 300
!...this check because getrel accepts + and - and no sign is allowed
   if(.not.(ch1.ge.'0' .and. ch1.le.'9') .and. ch1.ne.'.') then
      write(*,*)'ct1xfn 66:',ip,ch1,'>',string(ip-5:ip+5),'<'
      gx%bmperr=4001
      goto 1000
   endif
!   write(*,202)ip,string(1:ip+5)
202 format('Expect real: ',i5,' >',a,'< ')
   call getrel(string,ip,val)
   if(buperr.ne.0) then
      gx%bmperr=buperr
      goto 1000
   endif
! looking for 0*fun bug
   if(val.eq.zero) zeroc=.TRUE.
!   write(*,*)'ct1xfn 1: ',nterm,ip,val
!...if nterm>0 and coeff(nterm)=0 then store this coefficent there
   if(nterm.gt.0 .and. coeff(nterm).eq.zero) then
      coeff(nterm)=sign*val
!      write(*,*)'ct1xfn 2: ',nterm,val,coeff(nterm)
   else
      nterm=nterm+1
      if(nterm.gt.mterm) then
         gx%bmperr=4000; goto 1000
      endif
      coeff(nterm)=sign*val
      sign=one
!      write(*,*)'ct1xfn 3: ',nterm,val,coeff(nterm)
      do i=1,5
         koder(i,nterm)=0
      enddo
   endif
   goto 400
!
!...unsigned symbol, first character at ip
300 continue
   symbol=' '
   call ct1getsym(string,ip,symbol)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,*)'ct1xfn 5: ',nterm,ip,coeff(nterm)
!...one can have a symbol as first part, then create a term
!     otherwise symbols are usually part of a term already created
   if(nterm.eq.0) then
      nterm=1
      coeff(nterm)=one
      do i=1,5
         koder(i,nterm)=0
      enddo
   elseif(coeff(nterm).eq.zero) then
! this can happen if one has no coefficient in front of a function!!!
!      write(*,*)'ct1xfn 5A: ',nterm,ip,coeff(nterm)
      if(.not.zeroc) coeff(nterm)=sign*one
!      write(*,*)'ct1xfn 5B: ',nterm,ip,coeff(nterm)
   endif
!...check if T or P
   if(symbol(1:2).eq.tsym) then
      if(string(ip:ip).eq.'^' .or. string(ip:ip+1).eq.'**') then
         ip=ip+1
         if(string(ip:ip).eq.'*') ip=ip+1
         call ct1power(string,ip,ipower)
         if(gx%bmperr.ne.0) goto 1000
      else
         ipower=1
      endif
      koder(1,nterm)=koder(1,nterm)+ipower
      goto 400
   elseif(symbol(1:2).eq.psym) then
! allow powers as ^ or **
      if(string(ip:ip).eq.'^' .or. string(ip:ip+1).eq.'**') then
         ip=ip+1
         if(string(ip:ip).eq.'*') ip=ip+1
         call ct1power(string,ip,ipower)
         if(gx%bmperr.ne.0) goto 1000
      else
         ipower=1
      endif
      koder(2,nterm)=koder(2,nterm)+ipower
      goto 400
   endif
!...check if unary operator
   do j=1,nunary
      if(symbol(1:6).eq.unary(j)) goto 380
   enddo
! here search tpfuns for symbols, there are freetpfun-1 of them
   do jss=1,freetpfun-1
      if(symbol.eq.tpfuns(jss)%symbol) goto 350
   enddo
!...unknown new symbol
   if(fromtdb) then
! if we are reading a TDB file allow references to unknown functions
! We will scan for un-entered TPfuns later
!      write(*,*)'Unknown symbol to be entered later: ',symbol
      call enter_tpfun_dummy(symbol)
   else
! otherwise give error message
!      write(*,*)'TPFUN Unknown symbol: ',symbol,freetpfun-1
      gx%bmperr=4002; goto 1000
   endif
! we have found the symbol
350 continue
   if(koder(5,nterm).ne.0) then
! two symbols multipled with each other
      if(koder(3,nterm).ne.0) then
         write(*,*)'too many symbols in one term: ',koder(3,nterm)
         gx%bmperr=4022; goto 1000
      else
! set new function in koder(3,nterm), otherwise written in oposite order
         koder(3,nterm)=1000+jss
      endif
   else
      koder(5,nterm)=jss
   endif
   goto 400
!...unary function must be follwed by (
380 continue
   ch1=string(ip:ip)
   if(ch1.ne.'(') then
      gx%bmperr=4003
      goto 1000
   else
      ip=ip+1
      levelp=levelp+1
      koder(4,nterm)=levelp
      if(koder(5,nterm).ne.0) then
! this is like R*T*LN(1E-5*P), save link to R in koder(3,nterm)
         if(koder(3,nterm).ne.0) then
            write(*,*)'too many symbols in one term: ',koder(3,nterm)
            gx%bmperr=4022; goto 1000
         elseif(koder(5,nterm).lt.0) then
            write(*,*)'two unary functions in one term: ',koder(3,nterm)
            gx%bmperr=4023; goto 1000
         else
            koder(3,nterm)=1000+koder(5,nterm)
         endif
      endif
      koder(5,nterm)=-j
!...new term for argument of unary function, set coefficint to zero
!     to mark that none has been found.
      nterm=nterm+1
      if(nterm.gt.mterm) then
         write(*,*)'ct1xfn 8: ',nterm,mterm,ip
         gx%bmperr=4000
         goto 1000
      endif
      coeff(nterm)=zero
      do i=1,5
         koder(i,nterm)=0
      enddo
      goto 100
   endif
!
!...after a factor of a term: ),operator *, ^, +, - (division / not allowed)
400 continue
   if(eolch(string,ip)) goto 800
   ch1=string(ip:ip)
!...+ or - means new term
   if(ch1.eq.'-' .or. ch1.eq.'+') goto 100
   sign=one
   after: if(ch1.eq.')') then
      koder(4,nterm)=levelp
      if(levelp.eq.0) then
         gx%bmperr=4004
         goto 1000
      endif
      levelp=levelp-1
      ip=ip+1
      goto 400
   elseif(ch1.eq.'^') then
      ip=ip+1
      call ct1power(string,ip,ipower)
      if(koder(3,nterm).ne.0) then
! several symbols or unary and power, too complicated
         gx%bmperr=4024; goto 1000
      endif
      koder(3,nterm)=ipower
      goto 400
   elseif(ch1.eq.'*') then
      ip=ip+1
      ch1=string(ip:ip)
      if(ch1.eq.'*') then
         ip=ip+1
         call ct1power(string,ip,ipower)
         if(koder(3,nterm).ne.0) then
! several symbols or unariy and power, too complicated
            gx%bmperr=4024; goto 1000
         endif
         koder(3,nterm)=ipower
         goto 400
      else
! multiply symbol with something ....
!         write(*,*)'ct1xfn 4: ',nterm,ip,coeff(nterm)
654      format(a,i5,'"',a,'"',i5)
      endif
      goto 200
   elseif(ch1.eq.';') then
      goto 900
   endif after
     write(*,777)'error ct1xfn: ',ch1,ip,string(ip-10:ip+10)
777 format(a,'>',a,'<',i3,'>',a,'<')
   gx%bmperr=4005
   goto 1000
! no more characters, check!!
800 continue
!
!...; or no more characters, expression finished, check!!
900 continue
   if(levelp.gt.0) then
      gx%bmperr=4006
      goto 1000
   endif
   nc=nterm
990 format('ct1xfn 99> ',1PE15.6,5I7)
1000 continue
   return
 end subroutine ct1xfn

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine ct1getsym(string,ip,symbol)
!...extracts an symbol
!   implicit double precision (a-h,o-z)
   implicit none
   integer ip
   character string*(*),symbol*(*)
!\end{verbatim} %+
   integer, parameter :: lenfnsym=16
   integer jp
   character ch1*1,chs*1,localsym*(lenfnsym)
! these 2 functions are declared in METLIB and no type decration needed here
   jp=0
   localsym=' '
   symbol=' '
100 continue
   ch1=biglet(string(ip:ip))
!      write(6,*)'ct1getsym 2 >',ch1,'<',ip
   if((ch1.ge.'A' .and. ch1.le.'Z') .or. &
        (jp.gt.0 .and. ch1.eq.'_') .or. &
        (jp.gt.0 .and. (ch1.ge.'0' .and. ch1.le.'9'))) then
      jp=jp+1
! ignore characters after length of localsym
      if(jp.le.len(localsym)) then
         localsym(jp:jp)=ch1
      endif
      ip=ip+1
      goto 100
   endif
   symbol=localsym
1000 return
 end subroutine ct1getsym

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine ct1power(string,ip,ipower)
!...extracts an integer power possibly surrounded by ( )
   implicit none
   integer ip,ipower
   character string*(*)
!\end{verbatim} %+
   integer ich,isig,lp,jp
   character ch1*1,chs*1
   lp=0
   isig=1
   ipower=0
100 continue
   ch1=string(ip:ip)
   if(ch1.eq.'(') then
      if(lp.gt.0) then
         gx%bmperr=4007
         goto 1000
      elseif(ipower.ne.0) then
         gx%bmperr=4008
         goto 1000
      endif
      jp=ip+1
      if(eolch(string,jp)) then
         gx%bmperr=4009
         goto 1000
      endif
      chs=string(jp:jp)
      if(chs.eq.'-') then
!...mark ( and save sign, update ip (incremented below)
         lp=1
         isig=-1
         ip=jp
      endif
   elseif(ch1.eq.')') then
      if(ipower.ne.0) then
!...the ) can belong to other parts of the expression
         if(lp.eq.1) then
            ip=ip+1
            lp=0
         endif
         goto 900
      endif
      gx%bmperr=4010
      goto 1000
   elseif(ch1.ge.'0' .and. ch1.le.'9') then
      ich=ichar(ch1)-ichar('0')
      ipower=10*ipower+ich
   else
      goto 900
   endif
   ip=ip+1
   if(ipower.gt.100) then
      gx%bmperr=4011
      goto 1000
   endif
   goto 100
!
900 if(lp.gt.0) then
      gx%bmperr=4012
      goto 1000
   endif
   ipower=isig*ipower
1000 return
 end subroutine ct1power

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine ct1mfn(symbol,nranges,tlimits,lokexpr,lrot)
!...creates a root record with name symbol and temperature ranges
! highest T limit is in tlimits(nranges+1)
!   implicit double precision (a-h,o-z)
   implicit none
   integer nranges,lrot
   character*(*) symbol
   TYPE(tpfun_expression), dimension(*) :: lokexpr
   real tlimits(*)
!\end{verbatim} %+
   integer ir
   character name*16
   lrot=freetpfun
!   write(*,*)'ct1mfn: ',freetpfun
!   write(*,*)'ct1mfn: ',lrot,tpfuns(lrot)%nextfree
   if(lrot.gt.0) then
      freetpfun=tpfuns(lrot)%nextfree
      tpfuns(lrot)%nextfree=0
   else
! no more tpfun records
      write(*,*)'No more space for TP functions: ',size(tpfuns)
      gx%bmperr=4014; goto 1000
   endif
   allocate(tpfuns(lrot)%limits(nranges))
   allocate(tpfuns(lrot)%funlinks(nranges))
   do ir=1,nranges
      tpfuns(lrot)%limits(ir)=tlimits(ir)
! should this be an assignment or setting a link?
      tpfuns(lrot)%funlinks(ir)=lokexpr(ir)
   enddo
   tpfuns(lrot)%hightlimit=tlimits(nranges+1)
   tpfuns(lrot)%noofranges=nranges
! save name as upper case
   name=symbol
   call capson(name)
   tpfuns(lrot)%symbol=name
1000 continue
   return
 end subroutine ct1mfn

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine ct2mfn(symbol,nranges,tlimits,lokexpr,lrot)
!...stores a TPfun in an existing lrot record with name symbol
! and temperature ranges, highest T limit is in tlimits(nranges+1)
   implicit none
   integer nranges,lrot
   character*(*) symbol
   TYPE(tpfun_expression), dimension(*) :: lokexpr
   real tlimits(*)
!\end{verbatim} %+
   integer ir
   character name*16
   if(lrot.gt.0 .and. lrot.lt.freetpfun .and. &
        btest(tpfuns(lrot)%status,TPNOTENT)) then
      if(tpfuns(lrot)%noofranges.gt.0) then
         write(*,*)'This TPfun has already been entered ...',symbol
         gx%bmperr=7777; goto 1000
      endif
   else
! illegal value of lrot
      gx%bmperr=7777; goto 1000
   endif
   allocate(tpfuns(lrot)%limits(nranges))
   allocate(tpfuns(lrot)%funlinks(nranges))
   do ir=1,nranges
      tpfuns(lrot)%limits(ir)=tlimits(ir)
! should this be an assignment or setting a link?
      tpfuns(lrot)%funlinks(ir)=lokexpr(ir)
   enddo
   tpfuns(lrot)%hightlimit=tlimits(nranges+1)
   tpfuns(lrot)%noofranges=nranges
! clear the bit that this TPFUN is not entered
   tpfuns(lrot)%status=ibclr(tpfuns(lrot)%status,TPNOTENT)
!   write(*,*)'Clearing noentered bit: ',lrot,tpfuns(lrot)%symbol
! name already stored
!   name=symbol
!   call capson(name)
!   tpfuns(lrot)%symbol=name
1000 continue
   return
 end subroutine ct2mfn

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine ct1mexpr(nc,coeff,koder,lrot)
!...makes a datastructure of an expression. root is returned in lrot
!   implicit double precision (a-h,o-z)
   implicit none
   integer nc,koder(5,*)
   TYPE(tpfun_expression), pointer :: lrot
   TYPE(tpfun_expression), pointer :: noexpr
   double precision coeff(*)
!\end{verbatim} %+
   integer i
   if(nc.le.0) then
      nullify(lrot)
      goto 1000
   endif
! allocate an expression record and then allocate all arrays
   allocate(lrot)
   lrot%noofcoeffs=nc
   allocate(lrot%coeffs(nc))
   allocate(lrot%tpow(nc))
   allocate(lrot%ppow(nc))
   allocate(lrot%wpow(nc))
   allocate(lrot%plevel(nc))
   allocate(lrot%link(nc))
! store data
   save2: do i=1,nc
      lrot%coeffs(i)=coeff(i)
      lrot%tpow(i)=koder(1,i)
      lrot%ppow(i)=koder(2,i)
      lrot%wpow(i)=koder(3,i)
      lrot%plevel(i)=koder(4,i)
      lrot%link(i)=koder(5,i)
   enddo save2
1000  return
 end subroutine ct1mexpr

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine ct1efn(inrot,tpval,val,tpres)
!...evaluates a datastructure of an expression. Value returned in val
!     inrot is root expression tpfunction record
!     tpval is valuse of T and P, symval is values of symbols
! first and second derivatives of T and P also calculated and returned
! in order F, F.T, F.P, F.T.T, F.T.P, F.P.P
!
! if function already calculated one should never enter this subroutine
!
! It can call "itself" by reference to another TP function and for
! that case one must store results in levels.
   implicit none
   double precision val(6),tpval(*)
   TYPE(tpfun_expression), pointer :: inrot
   TYPE(tpfun_parres), dimension(:), pointer :: tpres
!\end{verbatim}
   integer mlev,level,jpow,link2
   double precision mini
   parameter (mlev=10,mini=1.0D-8)
   TYPE tpfun_nest
      TYPE(tpfun_nest), pointer :: previous
      TYPE(tpfun_expression), pointer :: exprot
      integer savenc,saveic,savelink4,level,savetp
      double precision saveval(6)
   end TYPE tpfun_nest
   TYPE(tpfun_nest), pointer :: temp,topsave
   TYPE(tpfun_expression), pointer :: exprot,nyrot
   double precision symval(6),sym,dsymdp,dsymdt
   double precision sym1,sym2,dsym1dt,dsym2dt,dsym1dp,dsym2dp
   double precision ff,dfdt,dfdp,d2fdt2,d2fdtdp,d2fdp2,cc
   double precision gg,dgdt,dgdp,d2gdt2,d2gdtdp,d2gdp2,cc1
   integer i,ic,mc,lrot,ipow,link,nc,unfun,itpow,tprot,link3,link4
   double precision t0,breakfun(6)
   integer becareful
!
   val=zero
   level=0
   link4=0
   exprot=>inrot
   nullify(topsave)
   becareful=0
!   write(*,*)'ct1efn 0: ',lrot,tpval(1)
!-----------------------------------------------
! return here for a linked function
100 continue
   if(.not.associated(exprot)) then
! this is not an error, just return zero or the value of a constant !!!
      goto 900
   endif
!------------------------------------------
   ic=0
   nc=exprot%noofcoeffs
!------------------------------------------------
!  return here for each new term and after evaluating linked symbols
200 continue
   eval: do while (ic.lt.nc)
      ic=ic+1
      cc=exprot%coeffs(ic)
      if(cc.eq.zero) cycle eval
      ipow=exprot%tpow(ic)
      if(ipow.ne.0) then
         ff=cc*tpval(1)**ipow
         dfdt=cc*ipow*tpval(1)**(ipow-1)
         dfdp=zero
         d2fdt2=cc*ipow*(ipow-1)*tpval(1)**(ipow-2)
         d2fdtdp=zero
         d2fdp2=zero
      else
         ff=cc
         dfdt=zero
         dfdp=zero
         d2fdt2=zero
         d2fdtdp=zero
         d2fdp2=zero
      endif
      ipow=exprot%ppow(ic)
      if(ipow.ne.0) then
! calculate backwards not to destroy value of ff
         d2fdp2=ff*ipow*(ipow-1)*tpval(2)**(ipow-2)
         d2fdtdp=dfdt*ipow*tpval(2)**(ipow-1)
         d2fdt2=d2fdt2*tpval(2)**ipow
         dfdp=ff*ipow*tpval(2)**(ipow-1)
         dfdt=dfdt*tpval(2)**ipow
         ff=ff*tpval(2)**ipow
      endif
!...power of symbols is handeled belown
!       ipow=exprot%wpow(ic)
      ipow=exprot%plevel(ic)
!...igore this at present, should never be set ...
!         if(ipow.ne.0) then
!           gx%bmperr=4017
!           goto 1000
!         endif
!>>>>>>>>>>>>> very uncertain code from here <<<<<<<<<<<<<<<<<<<<<<
      link=exprot%link(ic)
      link3=exprot%wpow(ic)
!       if(link.ne.0) write(*,201)'funev: ',lrot,ic,link,link3
201    format(a,10i5)
      if(link.lt.0 .and. link3.gt.1000) then
! if link is negative (unary funktion) and link3 is >1000 (link)
! we must evaluate link3 first
         link4=link3-1000
         if(abs(tpres(link4)%tpused(1)-tpval(1)).lt.&
              mini*tpres(link4)%tpused(1) .and. &
              abs(tpres(link4)%tpused(2)-tpval(2)).lt.&
              mini*tpres(link4)%tpused(2)) then
! function in link3-1000 is evaluated, multiply it with ff
            sym1=tpres(link4)%results(1)
            dsym1dt=tpres(link4)%results(2)
            dsym1dp=tpres(link4)%results(3)
            d2fdp2=sym1*d2fdp2+2.0D0*dsym1dp*dfdp+&
                 tpres(link4)%results(6)*ff
            d2fdtdp=sym1*d2fdtdp+dsym1dp*dfdt+dsym1dt*dfdp+&
                 tpres(link4)%results(5)*ff
            d2fdt2=sym1*d2fdt2+2.0D0*dsym1dt*dfdt+&
                 tpres(link4)%results(4)*ff
            dfdp=sym1*dfdp+dsym1dp*ff
            dfdt=sym1*dfdt+dsym1dt*ff
            ff=sym1*ff
!             write(*,202)'mulnk 0: ',lrot,ic,link,link4,0,sym1,ff
         else
! we must first evaluate link3, after that is done we come here again
! and take the else path below
            exprot%wpow(ic)=-1000+link
            link=link3-1000
            link4=link3
         endif
      elseif(link3.lt.-1000) then
! now we have evaluated the symbol, we must multiply that
! with ff and then evaluate the unary function
         d2fdp2=d2fdp2*symval(1)+2.0D0*dfdp*symval(3) &
              +ff*symval(6)
         d2fdtdp=d2fdtdp*symval(1)+dfdp*symval(2) &
              +dfdt*symval(3)+ff*symval(5)
         d2fdt2=d2fdt2*symval(1)+2.0D0*dfdt*symval(2) &
              +ff*symval(4)
         dfdp=dfdp*symval(1)+ff*symval(3)
         dfdt=dfdt*symval(1)+ff*symval(2)
         ff=ff*symval(1)
         exprot%wpow(ic)=link4
      endif
!-------------------------------------------------------------
      evlink: if(link.gt.0) then
! link to another symbol, extract its value and use chain rule
! extract the results from the symbol if already calculated
! if not calculated then do that and then recalculate this term
         linkif: if(abs(tpres(link)%tpused(1)-tpval(1)).lt.&
              mini*tpres(link)%tpused(1) .and. &
              abs(tpres(link)%tpused(2)-tpval(2)).lt.&
              mini*tpres(link)%tpused(2)) then
            jpow=exprot%wpow(ic)
!---------------------------------------------
            jpowif: if(jpow.gt.1000) then
! suck, two functions have to be multiplied ....
               link2=jpow-1000
               jpowev: if(abs(tpres(link2)%tpused(1)-tpval(1)).lt.&
                    mini*tpres(link2)%tpused(1) .and. &
                    abs(tpres(link2)%tpused(2)-tpval(2)).lt.&
                    mini*tpres(link2)%tpused(2)) then
! both functions are evaluated, multiply the two functions here
! one function is in tpres(link)%results, the other in tpres(link2)%results
                  sym1=tpres(link)%results(1)
                  dsym1dt=tpres(link)%results(2)
                  dsym1dp=tpres(link)%results(3)
                  sym2=tpres(link2)%results(1)
                  dsym2dt=tpres(link2)%results(2)
                  dsym2dp=tpres(link2)%results(3)
                  symval(6)=sym1*tpres(link2)%results(6)+&
                       2.0D0*dsym1dp*dsym2dp+&
                       tpres(link)%results(6)*sym2
                  symval(5)=sym1*tpres(link2)%results(5)+&
                       dsym1dp*dsym2dt+dsym1dt*dsym2dp+&
                       tpres(link)%results(5)*sym2
                  symval(4)=sym1*tpres(link2)%results(4)+&
                       2.0D0*dsym1dt*dsym2dt+&
                       tpres(link)%results(4)*sym2
                  symval(3)=sym1*dsym2dp+dsym1dp*sym2
                  symval(2)=sym1*dsym2dt+dsym1dt*sym2
                  symval(1)=sym1*sym2
               else
! function link2 must be evaluated, push and calculate
                  if(btest(tpfuns(link2)%status,TPCONST)) then
!                     write(*,*)'Link to a constant 1'
                     becareful=link2
                     nullify(nyrot)
                  else
                     call nested_tpfun(link2,tpval,nyrot)
                     if(gx%bmperr.ne.0) goto 1000
!                     write(*,*)'ct1efn nest 2: ',link2,nyrot
                  endif
! here we must push current values and start evaluating a new function nyrot
! when that has been done one must return here ... how??
! Well probably simplest by a new evaluation the same term again and when
! finding the link one takes the newly evaluated numbers !!!
! That means this function must save values in the tpfunction !!!
                  level=level+1
                  allocate(temp)
                  temp%previous=>topsave
                  topsave=>temp
                  nullify(temp)
                  topsave%exprot=>exprot
                  topsave%level=level
                  topsave%saveic=ic-1; topsave%savenc=nc
                  topsave%savetp=link2; topsave%savelink4=link4
                  topsave%saveval=val
                  val=zero
                  if(becareful.gt.0) then
! save the constant value in val(1), then jump to 900
                     val(1)=tpfuns(becareful)%limits(1)
                     becareful=0
                     goto 900
                  else
                     exprot=>nyrot
                     goto 100
                  endif
               endif jpowev
            elseif(jpow.eq.0 .or. jpow.lt.-1000) then
! jpow can be <-1000 if a symbol is multiplied with a unary function
! Here we just extract the values of the function
               do i=1,6
                  symval(i)=tpres(link)%results(i)
               enddo
            elseif(jpow.ne.0) then
! this symbol is raised to a power, use chain rule for derivatives backward
               sym=tpres(link)%results(1)
               dsymdt=tpres(link)%results(2)
               dsymdp=tpres(link)%results(3)
               symval(6)=jpow*(jpow-1)*sym**(jpow-2)*dsymdp**2+&
                    jpow*sym**(jpow-1)*tpres(link)%results(6)
               symval(5)=jpow*(jpow-1)*sym**(jpow-2)*dsymdp*dsymdt+&
                    jpow*sym**(jpow-1)*tpres(link)%results(5)
               symval(4)=jpow*(jpow-1)*sym**(jpow-2)*dsymdt**2+&
                    jpow*sym**(jpow-1)*tpres(link)%results(4)
               symval(3)=jpow*sym**(jpow-1)*dsymdp
               symval(2)=jpow*sym**(jpow-1)*dsymdt
               symval(1)=sym**jpow
            endif jpowif
         else
! one must evaluaste the function in link, it is recursive through eval_tpfun
! which will call ct1efn again but this is handelled automatically?????
! One should add some check that two TP functions does not call each other
! to infinite depth.  Same as done above
            if(btest(tpfuns(link)%status,TPCONST)) then
! the function is a constant!!
!               write(*,*)'Link to a constant 2'
               becareful=link
            else
               call nested_tpfun(link,tpval,nyrot)
               if(gx%bmperr.ne.0) goto 1000
            endif
! here we must push current values and start evaluating a new function nyrot
! when that has been done one must return here ... how??
! Well probably simplest: new evaluation of the same term again and when
! finding the link one takes the evaluated numbers !!!
! That means this function must save values in the tpfunction !!!
            level=level+1
            allocate(temp)
            temp%previous=>topsave
            topsave=>temp
            nullify(temp)
            topsave%exprot=>exprot
            topsave%level=level
            topsave%saveic=ic-1; topsave%savenc=nc
            topsave%savetp=link; topsave%savelink4=link4
            topsave%saveval=val
            val=zero
            if(becareful.gt.0) then
! save the constant value in val(1), then jump to 900
               val(1)=tpfuns(becareful)%limits(1)
               becareful=0
               goto 900
            else
               exprot=>nyrot
               goto 100
            endif
         endif linkif
! The symbol (or multiplied symbols) value in symval, apply chain rule
         d2fdp2=d2fdp2*symval(1)+2.0D0*dfdp*symval(3) &
              +ff*symval(6)
         d2fdtdp=d2fdtdp*symval(1)+dfdp*symval(2) &
              +dfdt*symval(3)+ff*symval(5)
         d2fdt2=d2fdt2*symval(1)+2.0D0*dfdt*symval(2) &
              +ff*symval(4)
         dfdp=dfdp*symval(1)+ff*symval(3)
         dfdt=dfdt*symval(1)+ff*symval(2)
         ff=ff*symval(1)
      elseif(link.lt.0) then
!------------------------------------------------------
! unary function, next term is argument, not very elegant ....
         unfun=link
         cc=exprot%coeffs(ic+1)
! cc should never be zero here, if so bug in the parser
         if(cc.eq.zero) then
            gx%bmperr=4018; goto 1000
         endif
         ipow=exprot%tpow(ic+1)
         if(ipow.ne.0) then
            gg=cc*tpval(1)**ipow
            dgdt=cc*ipow*tpval(1)**(ipow-1)
            dgdp=zero
            d2gdt2=cc*ipow*(ipow-1)*tpval(1)**(ipow-2)
            d2gdtdp=zero
            d2gdp2=zero
         else
            gg=cc
            dgdt=zero
            dgdp=zero
            d2gdt2=zero
            d2gdtdp=zero
            d2gdp2=zero
         endif
         ipow=exprot%ppow(ic+1)
         if(ipow.ne.0) then
            d2gdp2=gg*ipow*(ipow-1)*tpval(2)**(ipow-2)
            d2gdtdp=dgdt*ipow*tpval(2)**(ipow-1)
            d2gdt2=d2gdt2*tpval(2)**ipow
            dgdp=gg*ipow*tpval(2)**(ipow-1)
            dgdt=dgdt*tpval(2)**ipow
            gg=gg*tpval(2)**ipow
         endif
!...ignore these at present
         ipow=exprot%wpow(ic+1)
         if(ipow.ne.0) then
            gx%bmperr=4019
            goto 1000
         endif
         link2=exprot%link(ic+1)
         if(link2.gt.0) then
! link2, another symbol inside unary term, extract its value and use chain rule
! extract the results from the symbol if already calculated
            if(abs(tpres(link2)%tpused(1)-tpval(1)).lt.&
                 mini*tpres(link2)%tpused(1) .and. &
                 abs(tpres(link2)%tpused(2)-tpval(2)).lt.&
                 mini*tpres(link2)%tpused(2)) then
               symval=tpres(link2)%results
            else
! one must evaluaste another function, it is recursive through eval_tpfun
! which will call ct1efn again but this is handelled automatically?????
! One should add some check that two TP functions does not call each other
! to infinite depth
               if(btest(tpfuns(link2)%status,TPCONST)) then
! the function is a constant!!
!                  write(*,*)'The link is a constant 3'
                  becareful=link2
               else
                  call nested_tpfun(link2,tpval,nyrot)
                  if(gx%bmperr.ne.0) goto 1000
!                  write(*,*)'ct1efn nest 1: ',link2
               endif
! here we must push current values and start evaluating a new function nyrot
! when that has been done one must return here ... how??
! Well probably simplest ny evaluation the same term again and when
! finding the link one takes the evaluated numbers !!!
! That means this function must save values in the tpfunction !!!
               level=level+1
               allocate(temp)
               temp%previous=>topsave
               topsave=>temp
               nullify(temp)
               topsave%exprot=>exprot
               topsave%level=level
               topsave%saveic=ic-1; topsave%savenc=nc
               topsave%savetp=link2; topsave%savelink4=link4
               topsave%saveval=val
               val=zero
               if(becareful.gt.0) then
! save the constant value in val(1), then jump to 900
                  val(1)=tpfuns(becareful)%limits(1)
                  becareful=0
                  goto 900
               else
                  exprot=>nyrot
                  goto 100
               endif
            endif
            if(exprot%wpow(ic+1).ne.0) then
! it is illegal to have two symbols inside unary or power of symbol
               gx%bmperr=4016; goto 1000
            endif
! the value of the another symbol in symval.  use the chain rule
            d2gdp2=d2gdp2*symval(1)+2.0D0*dgdp*symval(3) &
                 +gg*symval(6)
            d2gdtdp=d2gdtdp*symval(1)+dgdp*symval(2) &
                 +dgdt*symval(3)+gg*symval(5)
            d2gdt2=d2gdt2*symval(1)+2.0D0*dgdt*symval(2) &
                 +gg*symval(4)
            dgdp=dgdp*symval(1)+gg*symval(3)
            dgdt=dgdt*symval(1)+gg*symval(2)
            gg=gg*symval(1)
         endif
! now combine term1 and term2 using chain rule. link values are
! -1: LOG,   -2: LN,    -3: EXP, -4: ERF, only LN and EXP implemented below
! -5: ABOVE, -6: BELOW are special functions for breakpoints
         evunfun: if(unfun.eq.-1) then
! ff=ff*Log10(gg) added by Sheng Yen Li
            if(gg.le.zero) then
               gx%bmperr=4020
               goto 1000
            endif
            d2fdp2=d2fdp2*log10(gg)+2.0d0*dfdp*dgdp/(gg*log(10d0)) &
                 -(ff*(dgdp/gg)**2)/log(10d0)+ff*d2gdp2/(gg*log(10d0))
            d2fdtdp=d2fdtdp*log10(gg)+dfdt*dgdp/(gg*log(10d0)) &
                 +dfdp*dgdt/(gg*log(10d0))-ff*dgdt*dgdp/((gg**2)*log(10d0)) &
                 +ff*d2gdtdp/(gg*log(10d0))
            d2fdt2=d2fdt2*log10(gg)+2.0d0*dfdt*dgdt/(gg*log(10d0)) &
                 -(ff*(dgdt/gg)**2)/log(10d0)+ff*d2gdt2/(gg*log(10d0))
            dfdp=dfdp*log10(gg)+ff*dgdp/(gg*log(10d0))
            dfdt=dfdt*log10(gg)+ff*dgdt/(gg*log(10d0))
            ff=ff*log10(gg)
         elseif(unfun.eq.-2) then
! ff=ff*LN(gg)
            if(gg.le.zero) then
               gx%bmperr=4020
               goto 1000
            endif
            d2fdp2=d2fdp2*log(gg)+2.0d0*dfdp*dgdp/gg &
                 -ff*(dgdp/gg)**2+ff/gg*d2gdp2
            d2fdtdp=d2fdtdp*log(gg)+dfdt*dgdp/gg+dfdp*dgdt/gg &
                 -ff*dgdt*dgdp/gg**2+ff*d2gdtdp/gg
            d2fdt2=d2fdt2*log(gg)+2.0d0*dfdt*dgdt/gg &
                 -ff*(dgdt/gg)**2+ff/gg*d2gdt2
            dfdp=dfdp*log(gg)+ff*dgdp/gg
            dfdt=dfdt*log(gg)+ff*dgdt/gg
            ff=ff*log(gg)
         elseif(unfun.eq.-3) then
! ff=ff*exp(gg)
            d2fdp2=exp(gg)*(d2fdp2+2.0D0*dfdp*dgdp+ff*d2gdp2+ff*(dgdp)**2)
            d2fdtdp=exp(gg)*(d2fdtdp+dfdt*dgdp+dfdp*dgdt+ff*d2gdtdp &
                 +ff*dgdt*dgdp)
            d2fdt2=exp(gg)*(d2fdt2+2.0D0*dfdt*dgdt+ff*d2gdt2+ff*(dgdt)**2)
            dfdp=exp(gg)*(dfdp+ff*dgdp)
            dfdt=exp(gg)*(dfdt+ff*dgdt)
            ff=ff*exp(gg)
         elseif(unfun.eq.-5) then
! above(T0): ff=ff*above(gg).  gg is t0, breakfun(1..6) are fun and derivatives
! NO P-derivatives means breakfun(3)=breakfun(5)=breakfun(6)=0
            call above_t0_calc(gg,tpval,breakfun)
            d2fdp2=breakfun(1)*d2fdp2
            d2fdtdp=breakfun(1)*d2fdtdp+dfdp*breakfun(2)
            d2fdt2=breakfun(1)*d2fdt2+ff*breakfun(4)+2.0D0*dfdt*breakfun(2)
            dfdp=breakfun(1)*dfdp
            dfdt=breakfun(1)*dfdt+ff*breakfun(2)
            ff=ff*breakfun(1)
         elseif(unfun.eq.-6) then
! below T0,  ff=ff*above(gg).  gg is t0, breakfun(1..6) are fun and derivatives
! NO P-derivatives means breakfun(3)=breakfun(5)=breakfun(6)=0
            call below_t0_calc(gg,tpval,breakfun)
            d2fdp2=breakfun(1)*d2fdp2
            d2fdtdp=breakfun(1)*d2fdtdp+dfdp*breakfun(2)
            d2fdt2=breakfun(1)*d2fdt2+ff*breakfun(4)+2.0D0*dfdt*breakfun(2)
            dfdp=breakfun(1)*dfdp
            dfdt=breakfun(1)*dfdt+ff*breakfun(2)
            ff=ff*breakfun(1)
         else
            gx%bmperr=4021
            goto 1000
         endif evunfun
         ic=ic+1
!----------------------------- end two-term unary function
!>>>>>>>>>>>>> very uncertain code above here <<<<<<<<<<<<<<<<<<<<<<
      else
! link=0, just continue
         continue
      endif evlink
! adding terms together
      val(1)=val(1)+ff
      val(2)=val(2)+dfdt
      val(3)=val(3)+dfdp
      val(4)=val(4)+d2fdt2
      val(5)=val(5)+d2fdtdp
      val(6)=val(6)+d2fdp2
   enddo eval
900 continue
! If level>1 save the values of the TP function.  The link to the
! address of TP function is in savetp(level)
   if(level.gt.0) then
! save calculated TP function values
      tprot=topsave%savetp
      do i=1,6
         tpres(tprot)%results(i)=val(i)
      enddo
      tpres(tprot)%tpused(1)=tpval(1)
      tpres(tprot)%tpused(2)=tpval(2)
! then unpack saved values of val and derivatives
      symval=val
      val=topsave%saveval
! POP the coefficients and the rest
      ic=topsave%saveic; nc=topsave%savenc
      link2=topsave%savetp; link4=topsave%savelink4
      exprot=>topsave%exprot
      topsave=>topsave%previous
      level=level-1
! restart from coefficient ic
       goto 200
    endif
!
1000 continue
   return
 end subroutine ct1efn

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine ct1wfn(exprot,tps,string,ip)
!...writes an expression into string starting at ip
!     lrot is an index to an tpexpr record
!   implicit double precision (a-h,o-z)
   implicit none
   character tps(2)*(*)
   character string*(*)
!\end{verbatim} %+
   integer, parameter :: levl=5,nunary=6
   integer, parameter :: lenfnsym=16
   integer koder(5,levl),ip,jus,is,kk,kpow,level,lpar,mult,nc,nos,i,ic
   double precision coeff(levl)
   character ch1*1,cht*1,extsym*(lenfnsym),unary(nunary)*6
   TYPE(tpfun_expression), pointer :: exprot
   DATA unary/'LOG   ','LN    ','EXP   ','ERF   ','ABOVE ','BELOW '/
!
   if(.not.associated(exprot)) then
      string(ip:ip+2)='0; '
      ip=ip+2
      goto 1000
   endif
   nc=exprot%noofcoeffs
   ic=0
   level=1
   lpar=0
!   write(*,*)'in ct1wfn',nc
200 ic=ic+1
   if(ic.gt.nc) goto 1000
   coeff(level)=exprot%coeffs(ic)
! bug in the expression parser ... (fixed??)
   if(coeff(level).eq.zero .and. nc.eq.1) then
! error in parser that such a function can exist, lrot should be zero
      string(ip:ip+1)='0 '
      ip=ip+1
      goto 1000
   elseif(coeff(level).eq.zero) then
      goto 200
   endif
   koder(1,level)=exprot%tpow(ic)
   koder(2,level)=exprot%ppow(ic)
   koder(3,level)=exprot%wpow(ic)
   koder(4,level)=exprot%plevel(ic)
   koder(5,level)=exprot%link(ic)
!
71    format(A,I5,1PE15.6,5I5)
   is=koder(5,level)
!   write(*,202)'ct1wfn: ',ic,ip,is,koder(1,level),string(1:ip)
!202 format(a,4i4,a)
   symbol: if(is.ne.0) then
!...reference to symbol or unary function, write coefficient only if not one
      if(abs(coeff(level)).ne.one) then
         call wrinum(string,ip,12,6,coeff(level))
         string(ip:ip)='*'
         ip=ip+1
         nos=0
         if(coeff(level).eq.zero) write(*,*)'ctwwfn ',ip,string(1:ip)
      elseif(coeff(level).eq.one) then
         nos=1
      else
         nos=-1
      endif
230    continue
      unaryfun: if(is.lt.0) then
!...write the T or P power before the unary function
         if(nos.eq.1) then
            string(ip:ip)='+'
            ip=ip+1
         elseif(nos.eq.-1) then
            string(ip:ip)='-'
            ip=ip+1
         endif
! there can be a symbol link in koder(3,level)
         if(koder(3,level).gt.1000) then
            jus=koder(3,level)-1000
            kk=len_trim(tpfuns(jus)%symbol)
            string(ip:ip+kk-1)=tpfuns(jus)%symbol
            ip=ip+kk
            string(ip:ip)='*'
            ip=ip+1
         endif
         call ct1wpow(string,ip,tps(1),1,koder(1,level))
         call ct1wpow(string,ip,tps(2),1,koder(2,level))
         kk=len_trim(unary(-is))
         string(ip:)=unary(-is)(1:kk)//'('
         ip=ip+kk+1
         lpar=koder(4,level)
!         write(*,*)'lpar: ',string(1:ip),' ',lpar
      else
! an external symbol, possibly a sign and power
         if(nos.eq.1) then
            string(ip:ip)='+'
            ip=ip+1
         elseif(nos.eq.-1) then
            string(ip:ip)='-'
            ip=ip+1
         endif
         kk=len_trim(tpfuns(is)%symbol)
         string(ip:ip+kk-1)=tpfuns(is)%symbol
         ip=ip+kk
         kpow=koder(3,level)
         if(kpow.gt.1000) then
! this is a link to another symbol, two symbols multiplied
            jus=kpow-1000
            kk=len_trim(tpfuns(jus)%symbol)
            string(ip:ip+kk)='*'//tpfuns(jus)%symbol
            ip=ip+kk+1
         elseif(kpow.lt.0) then
            kpow=-kpow
            if(kpow.gt.9) then
! power must be less than 99!!!
               ch1=char(ichar('0')+mod(kpow,10))
               cht=char(ichar('0')+kpow/10)
               string(ip:ip+6)='**(-'//cht//ch1//')'
               ip=ip+7
            else
               ch1=char(ichar('0')+kpow)
               string(ip:ip+5)='**(-'//ch1//')'
               ip=ip+6
            endif
         elseif(kpow.gt.0) then
! power must be less than 99!!!
            if(kpow.gt.9) then
               ch1=char(ichar('0')+mod(kpow,10))
               cht=char(ichar('0')+kpow/10)
               string(ip:ip+3)='**'//cht//ch1
               ip=ip+4
            else
               ch1=char(ichar('0')+kpow)
               string(ip:ip+2)='**'//ch1
               ip=ip+3
            endif
         endif
!...write the T or P power after the symbol and possible power
         call ct1wpow(string,ip,tps(1),-1,koder(1,level))
         call ct1wpow(string,ip,tps(2),-1,koder(2,level))
! fixing missing ) after unary function of symbol like exp(s1)
!         write(*,*)'problem here??:',string(1:ip),' ',lpar
! We got one extra ) as lpar not reset below
         if(lpar.gt.0) then
            string(ip:ip)=')'
            ip=ip+1
            lpar=0
         endif
      endif unaryfun
      goto 200
   endif symbol
! no symbol or unary function, coefficient with possible powers
   if(coeff(level).ne.one) then
      call wrinum(string,ip,12,6,coeff(level))
      mult=-1
   else
! in the case of a single value exactly 1 without unary or T or P power
! the number was never written
!      write(*,203)'ct1wfn2: ',(koder(i,level),i=1,4),coeff(level)
203   format(a,4i4,1pe12.4)
      do i=1,4
         if(koder(i,level).ne.0) goto 219
      enddo
! without this the Inden magnetic function will miss its initial 1.0
      call wrinum(string,ip,2,0,coeff(level))
      goto 220
219   continue
! missing coefficient discovered by Mauro, as the coefficient is unity
! it is not written.  Check with -1 maybe sign problems?
!      call wrinum(string,ip,2,1,coeff(level))
      string(ip:ip)='+'
      ip=ip+1
220   continue
      mult=0
   endif
!...write the T or P power after the coefficient
   call ct1wpow(string,ip,tps(1),mult,koder(1,level))
   call ct1wpow(string,ip,tps(2),mult,koder(2,level))
   if(koder(4,level).eq.1) then
      string(ip:ip)=')'
      ip=ip+1
! lpar was not reset here causing an extra ) later in expression ...
      lpar=0
!      write(*,*)'lpar not reset?:',string(1:ip)
!      write(*,*)lpar,koder(4,level)
   endif
   goto 200
1000 return
 end subroutine ct1wfn

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine ct1wpow(string,ip,tps,mult,npow)
!...writes "ips" with a power if needed and a * before or after
!   implicit double precision (a-h,o-z)
   implicit none
   integer ip,mult,npow
   character string*(*),tps*(*)
!\end{verbatim}
   integer lentps
   if(npow.eq.0) goto 1000
   if(mult.lt.0) then
      string(ip:ip)='*'
      ip=ip+1
   endif
   lentps=len_trim(tps)
   string(ip:ip+lentps-1)=tps
   ip=ip+lentps
   if(npow.gt.9) then
      write(string(ip:ip+3),110)npow
      ip=ip+4
   elseif(npow.gt.1) then
      write(string(ip:ip+4),120)npow
      ip=ip+3
   elseif(npow.lt.-9) then
      write(string(ip:ip+6),140)npow
      ip=ip+7
   elseif(npow.lt.0) then
      write(string(ip:ip+5),150)npow
      ip=ip+6
   endif
   if(mult.gt.0) then
      string(ip:ip)='*'
      ip=ip+1
   endif
110 format('**',i2)
120 format('**',i1)
140 format('**(',i3,')')
150 format('**(',i2,')')
1000 return
 end subroutine ct1wpow

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine enter_tpfun_interactivly(cline,ip,longline,jp)
! interactive input of a TP expression
!   implicit double precision (a-h,o-z)
   implicit none
   integer ip,jp
   character cline*(*),longline*(*)
!\end{verbatim}
   character line*80,ch1*1
   integer nexpr,lsc,kkp
   double precision xx
   call gparrd('Low temperature limit: ',cline,ip,xx,2.9815D2,nohelp)
   if(buperr.ne.0) then
      buperr=0; longline=' 298.15 '
      jp=8
   else
      longline=' '
      jp=1
      call wrinum(longline,jp,8,0,xx)
      if(buperr.ne.0) goto 1000
      jp=jp+1
   endif
   nexpr=1
   lsc=1
!-----------------------------------------------
! return here for new expression in another range
115 continue
   call gparc('Give expression, end with ";":',cline,ip,6,line,';',nohelp)
   if(buperr.ne.0) then
      buperr=0; line=';'
   endif
120 continue
   longline(jp:)=line
   jp=len_trim(longline)+1
!   write(*,*)'tpfun: ',longline(1:jp)
! lsc is position after the ";" in any previous range
   if(index(longline(lsc:),';').le.0) then
      call gparc('&',cline,ip,6,line,';',nohelp)
      if(buperr.ne.0) then
         buperr=0; line=';'
      endif
      goto 120
   endif
150 continue
! make sure there is a ; at the end of each expression
!   write(*,*)'tpfun add ;'
   kkp=index(longline(nexpr:),';')
   if(kkp.le.0) then
      kkp=len_trim(longline)
      longline(kkp+1:)='; '
      jp=kkp+3
      nexpr=jp
   else
      nexpr=kkp+1
   endif
! lsc is position of ; for previous range
   lsc=nexpr
   call gparrd('Upper temperature limit ',cline,ip,xx,6.0D3,nohelp)
   if(buperr.ne.0) then
      buperr=0; xx=6.0D3
   endif
   call wrinum(longline,jp,8,0,xx)
   if(buperr.ne.0) goto 1000
   call gparcd('Any more ranges',cline,ip,1,ch1,'N',nohelp)
   if(ch1.eq.'n' .or. ch1.eq.'N') then
      longline(jp:)=' N'
      jp=jp+3
   else
      longline(jp:)=' Y'
      jp=jp+3
      goto 115
   endif
1000 continue
   return
 end subroutine enter_tpfun_interactivly

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!-\begin{verbatim}
! subroutine tpfun_deallocate
! deallocates all arrays associated with a TP function
!-\end{verbatim}
!   implicit double precision (a-h,o-z)
!   implicit none
!   integer j,nr
!   do j=1,freetpfun-1
!      nr=tpfuns(j)%noofranges
!      if(nr.gt.0) then
!         deallocate(tpfuns(j)%funlinks)
!         deallocate(tpfuns(j)%limits)
!         deallocate(tpfuns(j)%funlinks)
!      endif
!   enddo
!1000 continue
!   return
! end subroutine tpfun_deallocate

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine enter_tpfun_dummy(symbol)
! creates a dummy entry for a TP function called symbol, used when entering 
! TPfuns from a TDB file where they are not in order
   implicit none
   character*(*) symbol
!\end{verbatim}
! set the TPNOTENT bit of this symbol
   integer lrot
   character name*16
   lrot=freetpfun
   if(lrot.gt.0) then
      freetpfun=tpfuns(lrot)%nextfree
      tpfuns(lrot)%nextfree=0
   else
      write(*,*)'No space for TP functions: ',size(tpfuns)
      gx%bmperr=4014; goto 1000
   endif
   tpfuns(lrot)%noofranges=0
   name=symbol
   call capson(name)
   tpfuns(lrot)%symbol=name
   tpfuns(lrot)%status=ibset(tpfuns(lrot)%status,TPNOTENT)
1000 continue
   return
 end subroutine enter_tpfun_dummy

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine enter_tpfun(symbol,text,lrot,fromtdb)
! creates a data structure for a TP function called symbol with several ranges
! text is whole expression
! lrot is returned as index.  If fromtdb is FALSE and lrot<0 it is a new
!                             expression for an old symbol
! if fromtdb is TRUE references to unknown functions are allowed
! default low temperature limit is 298.16; high 6000
   implicit none
   integer lrot
   character*(*) text,symbol
   logical fromtdb
!\end{verbatim}
! max number of ranges, max number of coefficents in each range
   integer, parameter :: mrange=20,mc=15
   integer jss,nc,ip,nrange
   real tlim(mrange)
   double precision coeff(mc),val
   integer koder(5,mc)
   TYPE(tpfun_expression) :: links(mrange)
   TYPE(tpfun_expression), pointer :: ltpexpr
   character ch1*1,lsym*(lenfnsym)
   logical already
! check if function already entered, there are freetpfun-1 of them
! ignore functions that start with a "_" as they are parameters
!   lrot=0
   already=.FALSE.
   if(symbol(1:1).ne.'_') then
      lsym=symbol
      call capson(lsym)
      do jss=1,freetpfun-1
!         write(*,17)jss,lsym,tpfuns(jss)%symbol
17       format('enter_tpfun: ',i5,' >,'a,'=',a,'?')
         if(lsym.eq.tpfuns(jss)%symbol) then
            if(btest(tpfuns(jss)%status,TPNOTENT)) then
! function name already eneterd, now enter expression, this is from TDB files
               lrot=jss; already=.TRUE.; goto 18
            else
!               write(*,*)'amend tpfun: ',fromtdb,lrot
               if(.NOT.fromtdb .and. lrot.lt.0) then
! this is an AMEND TPFUN, delete old expression to be able to store a new
                  lrot=jss; already=.TRUE.
                  nrange=tpfuns(lrot)%noofranges
!                  write(*,*)'Deallocating: ',lrot,nrange
                  deallocate(tpfuns(lrot)%limits)
                  deallocate(tpfuns(lrot)%funlinks)
                  tpfuns(lrot)%noofranges=0
                  tpfuns(lrot)%status=ibset(tpfuns(lrot)%status,TPNOTENT)
! we have to clear the stored values! But those are stored in equil.rec
                  nrange=0; goto 18
               else
                  gx%bmperr=4026; goto 1000
               endif
            endif
         endif
      enddo
   endif
!
   lrot=0
18 continue
! low T limit
   ip=1
   call getrel(text,ip,val)
   if(buperr.ne.0) then
! A , has been used to select default low temperature limit
      write(*,*)'Warning low temperature limit set to 298.15 K'
      val=298.15; buperr=0
      write(*,19)ip,text(1:72)
19    format('TPFUN: ',i3,' >',a)
   endif
   tlim(1)=val
   nrange=0
   ch1='Y'
! parse and store expression for each temperature range
   ranges: do while(ch1.eq.'Y')
      nrange=nrange+1
      if(nrange.gt.mrange) then
         gx%bmperr=4025; goto 1000
      endif
      nc=mc
      call ct1xfn(text,ip,nc,coeff,koder,fromtdb)
      if(gx%bmperr.ne.0) goto 1000
      call ct1mexpr(nc,coeff,koder,ltpexpr)
      links(nrange)=ltpexpr
! bypass final ; of expression
      ip=ip+1
      call getrel(text,ip,val)
      if(buperr.ne.0) then
! acceppt a , for default ...
         if(text(ip:ip).eq.',') then
            val=6.0D3; buperr=0
         else
            write(*,27)buperr,ip,text(1:ip+5)
27          format(' *** Error in enter_tpfun 2: ',i5,', position ',i5/&
                 '>',a,'<')
         endif
      endif
      tlim(nrange+1)=val
      if(.not.eolch(text,ip)) then
         ch1=biglet(text(ip:ip))
         ip=ip+1
      endif
   enddo ranges
   if(already) then
! a function symbol already entered, lrot is location
      call ct2mfn(symbol,nrange,tlim,links,lrot)
   else
! a new function record will be allocated
      call ct1mfn(symbol,nrange,tlim,links,lrot)
   endif
1000 continue
   return
 end subroutine enter_tpfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine nested_tpfun(lrot,tpval,nyrot)
! called from ct1efn when a it calls another TP function that must be
! evaluated.  nyrot is the link to the ct1efn in the correct range
!   implicit double precision (a-h,o-z)
   implicit none
   integer lrot
   double precision tpval(2)
   TYPE(tpfun_expression), pointer :: nyrot
! use lowest range for all T values lower than first upper limit
! and highest range for all T values higher than the next highest limit
! one should signal if T is lower than lowest limit or higher than highest
! used  saved reults if same T and P
!\end{verbatim} %+
   integer nr,ns
   nullify(nyrot)
   if(lrot.le.0) goto 1000
   nr=tpfuns(lrot)%noofranges
   if(nr.eq.0) then
! this is the case for constants! Does this work??
      if(btest(tpfuns(lrot)%status,TPCONST)) then
         write(*,*)'nested constant: ',nr,lrot
      else
         write(*,*)'A never never error evaluation a TP function',lrot
         write(*,*)'Function name: ',tpfuns(lrot)%symbol
         gx%bmperr=6666; goto 1000
      endif
   elseif(nr.eq.1) then
      nyrot=>tpfuns(lrot)%funlinks(1)
   else
      ns=1
      do while(ns.lt.nr)
         if(tpval(1).lt.tpfuns(lrot)%limits(ns+1)) then
            nyrot=>tpfuns(lrot)%funlinks(ns)
            goto 900
         endif
         ns=ns+1
      enddo
      nyrot=>tpfuns(lrot)%funlinks(nr)
   endif
900 continue
1000 continue
   return
 end subroutine nested_tpfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 logical function compare_abbrev(name1,name2)
! returns TRUE if name1 is an abbreviation of name2
! termintaes when a space is found in name1
! each part between _ or - can be abbreviated from the left
! case insensitive. Only 36 first characters compared
   implicit none
   character*(*) name1,name2   
!\end{verbatim} %+
   integer, parameter :: maxl=36
   integer jp,ip,noabbr
   character ch1*1
   character (len=maxl) :: lname1,lname2
   lname1=name1; lname2=name2
   call capson(lname1)
   call capson(lname2)
   compare_abbrev=.FALSE.
   noabbr=0
   jp=1
   bigloop: do ip=1,36
      ch1=lname1(ip:ip)
      if(ip.gt.1 .and. ch1.eq.' ') goto 900
      if(ch1.eq.'-') ch1='_'
      if(ch1.eq.lname2(jp:jp)) goto 300
      if(ch1.eq.'_' .or. ch1.eq.'-') then
200       continue
         if(jp.eq.maxl) goto 1000
         jp=jp+1
         if(lname2(jp:jp).eq.'_') goto 300
         if(lname2(jp:jp).eq.' ') goto 1000
         goto 200
      endif
      goto 1000
300    continue
      jp=jp+1
310    continue
   enddo bigloop
900 continue
   compare_abbrev=.TRUE.
1000 continue
   return
 end function compare_abbrev

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine below_t0_calc(t0,tpval,fun)
! calculates exp(20(1-t/t0))/(1+exp(20(1-t/t0))
! At t<<t0 K function is unity, at t>>t0 function is zero
   implicit none
   double precision t0,tpval(2),fun(6)
!\end{verbatim} %+
   double precision arg,expa
   double precision, parameter :: ffix=2.0D1
   if(t0.le.zero) then
      write(*,*)'temperature breakpoint below zero'
      goto 1000
   endif
   arg=ffix*(one-tpval(1)/t0)
   expa=exp(arg)
! F, F.T, T.P, F.T.T, F.T.P, F.P.P
!   fun(1)=-t0/ffix*log(1+expa)
!   fun(2)=expa/(one+expa)
!   fun(4)=-ffix*expa/(t0*(one+expa)**2)
   fun(1)=expa/(one+expa)
   fun(2)=-ffix*expa/(t0*(one+expa)**2)
   fun(4)=(ffix/t0)**3*expa*(one-expa)/(one+expa)**3
   fun(3)=zero
   fun(5)=zero
   fun(6)=zero
1000 continue
   return
 end subroutine below_t0_calc

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine above_t0_calc(t0,tpval,fun)
! calculates exp(20(1-t/t0))/(1+exp(20(1-t/t0))-1, at t>>t0 it is unity
   implicit none
   double precision t0,tpval(2),fun(6)
!\end{verbatim}
   call below_t0_calc(t0,tpval,fun)
   fun(1)=one-fun(1)
   fun(2)=-fun(2)
   fun(4)=-fun(4)
1000 continue
   return
 end subroutine above_t0_calc

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine enter_optvars(firstindex)
! enter variables for optimization A00-A99
   implicit none
   integer firstindex
!\end{verbatim} %+
   character symbol*(lenfnsym)
   integer jss,symix,lrot
   symbol='A00 '
! check if any TP fun with name A00 already entered
   do jss=1,freetpfun-1
      if(symbol.eq.tpfuns(jss)%symbol) then
         write(kou,*)'Optimizing symbols already entered'
         goto 1000
      endif
   enddo
   firstindex=freetpfun
   do jss=1,100
! create TPfun symbols with names A00 to A99 with value 0.0D0
      lrot=freetpfun
      if(lrot.eq.0) then
         gx%bmperr=4104; goto 1000
      else
         freetpfun=tpfuns(lrot)%nextfree
         tpfuns(lrot)%nextfree=0
      endif
      allocate(tpfuns(lrot)%limits(1))
      allocate(tpfuns(lrot)%funlinks(1))
      tpfuns(lrot)%symbol=symbol
      tpfuns(lrot)%limits(1)=zero
! mark this is a single value and can be optimized
      tpfuns(lrot)%status=ibset(tpfuns(lrot)%status,TPCONST)
      tpfuns(lrot)%status=ibset(tpfuns(lrot)%status,TPOPTCON)
! increment symbol
      symix=ichar(symbol(3:3))-ichar('0')
      symix=symix+1
      if(symix.eq.10) then
         symbol(3:3)='0'
         symbol(2:2)=char(ichar(symbol(2:2))+1)
      else
         symbol(3:3)=char(ichar(symbol(3:3))+1)
      endif
!      write(*,*)'Next symbol created: ',symbol(1:4),lrot
   enddo
1000 continue
   return
 end subroutine enter_optvars

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine find_tpsymbol(name,type,value)
! enter variables 
   implicit none
! type=0 if function, 1 if variable, 2 if optimizing variable
   integer type
   character name*(lenfnsym)
   double precision value
!\end{verbatim} %+
   integer jss,symix,lrot
   character symbol*(lenfnsym)
   symbol=name
   call capson(symbol)
! check if any TP fun with name symbol exists
   type=0
   do jss=1,freetpfun-1
      if(symbol.eq.tpfuns(jss)%symbol) then
! found symbol
         if(btest(tpfuns(jss)%status,TPCONST)) then
            value=tpfuns(jss)%limits(1)
            if(btest(tpfuns(jss)%status,TPOPTCON)) then
               type=2
            else
               type=1
            endif
         endif
         goto 200
      endif
   enddo
! no such symbol
   gx%bmperr=7777
   type=-1
200 continue
1000 continue
   return
 end subroutine find_tpsymbol

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine enter_tpconstant(symbol,value)
! enter variables 
   implicit none
   character symbol*(lenfnsym)
   double precision value
!\end{verbatim} %+
   integer jss,symix,lrot
! check if any TP fun with name symbol already entered
   do jss=1,freetpfun-1
      if(symbol.eq.tpfuns(jss)%symbol) then
! symbol already exist, just change value unless it is an optimizing coeff.
         if(btest(tpfuns(jss)%status,TPOPTCON)) then
            write(*,*)'Not allowed to change optimizing coefficents'
            goto 1000
         else
            lrot=jss
            goto 200
         endif
      endif
   enddo
! create TPfun symbols with name symbol and value value
   lrot=freetpfun
   if(lrot.eq.0) then
      gx%bmperr=4104; goto 1000
   else
      freetpfun=tpfuns(lrot)%nextfree
      tpfuns(lrot)%nextfree=0
   endif
   allocate(tpfuns(lrot)%limits(1))
   allocate(tpfuns(lrot)%funlinks(1))
   call capson(symbol)
   tpfuns(lrot)%symbol=symbol
! mark this is a single value
   tpfuns(lrot)%status=ibset(tpfuns(lrot)%status,TPCONST)
200 continue
   tpfuns(lrot)%limits(1)=value
1000 continue
   return
 end subroutine enter_tpconstant

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine change_optcoeff(lrot,value)
! change value of optimizing coefficient.  lrot is index
   implicit none
   integer lrot
   double precision value
!\end{verbatim} %+
   if(lrot.le.0 .or. lrot.ge.freetpfun-1) then
      write(*,*)'Attempt to change non-existing constant',lrot
      gx%bmperr=7777; goto 1000
   endif
   if(.not.btest(tpfuns(lrot)%status,TPOPTCON)) then
      write(*,*)'Attempt to change non-existing coefficent',lrot
      gx%bmperr=7777; goto 1000
   endif
   tpfuns(lrot)%limits(1)=value
1000 continue
   return
 end subroutine change_optcoeff

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine get_value_of_constant_name(symbol,lrot,value)
! get value (and index) of a TP constant.  lrot is index
   implicit none
   integer lrot
   character symbol*(*)
   double precision value
!\end{verbatim} %+
   write(*,*)'get_value_of_constant_name not implemented yet'
!   value=tpfuns(lrot)%limits(1)
1000 continue
   return
 end subroutine get_value_of_constant_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine get_value_of_constant_index(lrot,value)
! get value of a TP constant at known lrot
   implicit none
   integer lrot
   double precision value
!\end{verbatim} %+
   if(lrot.le.0 .or. lrot.gt.freetpfun-1) then
      write(kou,*)'Constant index outside limits',lrot
   else
! unifished: check if it is really a constant ...
      value=tpfuns(lrot)%limits(1)
   endif
1000 continue
   return
 end subroutine get_value_of_constant_index

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine get_all_opt_coeff(values)
! get values of all optimizing coefficients
   implicit none
   double precision values(*)
!\end{verbatim} %+
   write(*,*)'Not yet implemeneted'
1000 continue
   return
 end subroutine get_all_opt_coeff

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine delete_all_tpfuns
! delete all TPFUNs.  No error if some are already deleted ...   
! note: tpres is deallocated when deleting equilibrium record
!\end{verbatim}
   implicit none
   integer lrot,nrex
   TYPE(tpfun_expression), pointer :: expr
!   write(*,*)'In delete_all_tpfuns'
   deallocate(tpfuns)
   goto 1000
! code below skipped as it created a lot of memory errors ...
   if(tpfun_expression_version.ne.1 .or. tpfun_root_version.ne.1 .or. &
        tpfun_parres_version.ne.1) then
      write(*,*)'Data structure error when deleting tpfuns',&
           tpfun_expression_version,tpfun_root_version,&
           tpfun_parres_version
      gx%bmperr=7777; goto 1000
   endif
   funloop: do lrot=1,freetpfun-1
      write(*,*)'TP Deleting TP function: ',lrot
!      if(tpfuns(lrot)%noofranges.eq.0) cycle
      if(tpfuns(lrot)%noofranges.eq.0) goto 200
      write(*,*)'TP deleting ranges 1-',tpfuns(lrot)%noofranges
      range: do nrex=1,tpfuns(lrot)%noofranges
         expr=>tpfuns(lrot)%funlinks(nrex)
         if(associated(expr)) then
            deallocate(expr%coeffs)
            deallocate(expr%tpow)
            deallocate(expr%ppow)
            deallocate(expr%wpow)
            deallocate(expr%plevel)
            deallocate(expr%link)
            deallocate(expr)
         else
            write(*,*)'TP delete; no expression? ',lrot,nrex
         endif
      enddo range
200   continue
      write(*,*)'TP deleting limits ',size(tpfuns(lrot)%limits)
      deallocate(tpfuns(lrot)%limits)
      write(*,*)'TP deleting funlinks ',size(tpfuns(lrot)%funlinks)
      deallocate(tpfuns(lrot)%funlinks)
   enddo funloop
   write(*,*)'TP finally deleting roots'
   deallocate(tpfuns)
1000 continue
   return
 end subroutine delete_all_tpfuns

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!-\begin{verbatim}
! subroutine tpfunlim(i1,i2)
! used when saving TP funs to files in GTP as highexpr and freexpr are private
! and they are saved at several places in the file (but that can be changed)
!   implicit none
!   integer i1,i2
!-\end{verbatim}
!   i1=highexpr
!   i2=freexpr
!   return
! end subroutine tpfunlim

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine tpfunsave(lut,form)
! save tpfundata on a file, unfinished
!   implicit double precision (a-h,o-z)
   implicit none
   integer lut
   logical form
!\end{verbatim} %+
   integer i,j,kx,nc,nr
   TYPE(tpfun_expression), pointer :: exprot
! save the number of functions to read
   if(form) then
      write(lut,17)freetpfun
17    format('Number of TPFUNS: ',i7)
   else
      write(lut)tpfun_expression_version,tpfun_root_version,&
           tpfun_parres_version,freetpfun
   endif
!----------- tpfunction root list
   do j=1,freetpfun-1
      if(form) then
         write(lut,240)j,tpfuns(j)%symbol,tpfuns(j)%noofranges,&
              tpfuns(j)%nextfree
         nr=tpfuns(j)%noofranges
         if(nr.gt.0) then
            write(lut,241)(tpfuns(j)%limits(i),i=1,nr)
            do kx=1,nr
               exprot=>tpfuns(j)%funlinks(kx)
! the expressions
               nc=exprot%noofcoeffs
               if(nc.gt.0) then
                  write(lut,231)(exprot%coeffs(i),i=1,nc)
                  write(lut,235)(exprot%tpow(i),exprot%ppow(i),&
                       exprot%wpow(i),exprot%plevel(i),exprot%link(i),i=1,nc)
!                  write(lut,232)(exprot%tpow(i),i=1,nc)
!                  write(lut,232)(exprot%ppow(i),i=1,nc)
!                  write(lut,232)(exprot%wpow(i),i=1,nc)
!                  write(lut,232)(exprot%plevel(i),i=1,nc)
!                  write(lut,232)(exprot%link(i),i=1,nc)
               endif
230            format('TPEXPRESSION ',i5,'; ',i10)
231            format(('&',5(1PE15.7)))
232            format(('&',5i15))
235            format(('&',25i3))
            enddo
         endif
         write(lut,243)tpfuns(j)%hightlimit
240      format('TPFUNREC ',i5,': ',a,2i10)
241      format(('&',10F7.2))
242      format(('&',10i7))
243      format('&',F10.2)
      else
! unformatted
         write(lut)tpfuns(j)%symbol,tpfuns(j)%noofranges,tpfuns(j)%nextfree
         nr=tpfuns(j)%noofranges
         if(nr.gt.0) then
            write(lut)(tpfuns(j)%limits(i),i=1,nr)
            do kx=1,nr
               exprot=>tpfuns(j)%funlinks(kx)
               write(lut)exprot%noofcoeffs
               nc=exprot%noofcoeffs
               if(nc.gt.0) then
                  write(lut)(exprot%coeffs(i),i=1,nc)
                  write(lut)(exprot%tpow(i),i=1,nc)
                  write(lut)(exprot%ppow(i),i=1,nc)
                  write(lut)(exprot%wpow(i),i=1,nc)
                  write(lut)(exprot%plevel(i),i=1,nc)
                  write(lut)(exprot%link(i),i=1,nc)
               endif
            enddo
         else
            
         endif
         write(lut)tpfuns(j)%hightlimit
      endif
   enddo
!----------- tpfunction result list should not be save here
!   if(form) then
!      write(lut,250)(j,tpres(j)%tpused,tpres(j)%results,j=1,freetpfun-1)
!250    format(('TPRES ',i5,'; ',2(1PE15.6)/'&',6F13.5))
!   else
!      write(lut)(tpres(j)%tpused,tpres(j)%results,j=1,freetpfun-1)
!   endif
!
1000 continue
   return
 end subroutine tpfunsave

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine save1tpfun(lut,form,jfun)
! save one tpfun (a parameter) with index jfun on a file
!   implicit double precision (a-h,o-z)
   implicit none
   integer lut,jfun
   logical form
!\end{verbatim} %+
   integer nr,i,kx,nc
   TYPE(tpfun_expression), pointer :: exprot
   double precision dummy
! jfun can be zero meaning a parameter that is zero   
  if(form) then
     if(jfun.eq.0) goto 1000
      write(lut,240)jfun,tpfuns(jfun)%symbol,tpfuns(jfun)%noofranges,&
           tpfuns(jfun)%nextfree
      nr=tpfuns(jfun)%noofranges
      if(nr.gt.0) then
         write(lut,241)(tpfuns(jfun)%limits(i),i=1,nr)
         do kx=1,nr
            exprot=>tpfuns(jfun)%funlinks(kx)
! the expressions
            nc=exprot%noofcoeffs
            if(nc.gt.0) then
               write(lut,231)(exprot%coeffs(i),i=1,nc)
               write(lut,235)(exprot%tpow(i),exprot%ppow(i),&
                    exprot%wpow(i),exprot%plevel(i),exprot%link(i),i=1,nc)
!               write(lut,232)(exprot%tpow(i),i=1,nc)
!               write(lut,232)(exprot%ppow(i),i=1,nc)
!               write(lut,232)(exprot%wpow(i),i=1,nc)
!               write(lut,232)(exprot%plevel(i),i=1,nc)
!               write(lut,232)(exprot%link(i),i=1,nc)
            endif
230         format('TPEXPRESSION ',i5,'; ',i10)
231         format(('&',5(1PE15.7)))
232         format(('&',5i15))
235         format(('&'25i3))
         enddo
      endif
      write(lut,243)tpfuns(jfun)%hightlimit
240   format('TPFUNREC ',i5,': ',a,2i10)
241   format(('&',10F7.2))
242   format(('&',10i7))
243   format('&',F10.2)
   else
! unformatted
      if(jfun.eq.0) then
         write(lut)0,'ZERO            ',0,0
         write(lut)dummy
         write(*,*)'Writing a zero parameter'
         goto 1000
      endif
      write(lut)jfun,tpfuns(jfun)%symbol,tpfuns(jfun)%noofranges,&
           tpfuns(jfun)%nextfree
      nr=tpfuns(jfun)%noofranges
      if(nr.gt.0) then
         write(lut)(tpfuns(jfun)%limits(i),i=1,nr)
         do kx=1,nr
            exprot=>tpfuns(jfun)%funlinks(kx)
            write(lut)exprot%noofcoeffs
            nc=exprot%noofcoeffs
            if(nc.gt.0) then
               write(lut)(exprot%coeffs(i),i=1,nc)
               write(lut)(exprot%tpow(i),i=1,nc)
               write(lut)(exprot%ppow(i),i=1,nc)
               write(lut)(exprot%wpow(i),i=1,nc)
               write(lut)(exprot%plevel(i),i=1,nc)
               write(lut)(exprot%link(i),i=1,nc)
            endif
         enddo
      endif
      write(lut)tpfuns(jfun)%hightlimit
   endif
1000 continue
   return
 end subroutine save1tpfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine read1tpfun(lut,jfun)
! read one unformatted tpfun (a parameter) with index jfun
!   implicit double precision (a-h,o-z)
   implicit none
   integer lut,jfun
!\end{verbatim} %+
   integer i,i2,kx,nc,nr
   TYPE(tpfun_expression), pointer :: exprot
   character*16 symbol
   double precision dummy
! jfun can be sero meaning no link to a TPFUN
   read(lut)jfun,symbol,nr,i2
   if(jfun.gt.0) then
      tpfuns(jfun)%symbol=symbol
      tpfuns(jfun)%noofranges=nr
      tpfuns(jfun)%nextfree=i2
   else
      write(*,*)'reading a zero function: ',symbol
   endif
!   write(*,17)jfun,tpfuns(jfun)%symbol,tpfuns(jfun)%noofranges
!17 format('read1tpfun: ',i5,2x,a,i7)
   if(nr.gt.0) then
      allocate(tpfuns(jfun)%limits(nr))
      allocate(tpfuns(jfun)%funlinks(nr))
      read(lut)(tpfuns(jfun)%limits(i),i=1,nr)
! read the expressions
      do kx=1,nr
         exprot=>tpfuns(jfun)%funlinks(kx)
         read(lut)exprot%noofcoeffs
         nc=exprot%noofcoeffs
         if(nc.gt.0) then
            allocate(exprot%coeffs(nc))
            allocate(exprot%tpow(nc))
            allocate(exprot%ppow(nc))
            allocate(exprot%wpow(nc))
            allocate(exprot%plevel(nc))
            allocate(exprot%link(nc))
            read(lut)(exprot%coeffs(i),i=1,nc)
            read(lut)(exprot%tpow(i),i=1,nc)
            read(lut)(exprot%ppow(i),i=1,nc)
            read(lut)(exprot%wpow(i),i=1,nc)
            read(lut)(exprot%plevel(i),i=1,nc)
            read(lut)(exprot%link(i),i=1,nc)
         endif
      enddo
   endif
   if(jfun.gt.0) then
      read(lut)tpfuns(jfun)%hightlimit
   else
      read(lut)dummy
   endif
1000 continue
   return
 end subroutine read1tpfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine tpfunread(lin,skip)
! read tpfundata from a file, if skip TRUE ignore if function already entered
   implicit none
   integer lin
   logical skip
!\end{verbatim}
   integer i,j,kx,nr,nc,ranges,free,lrot
   character symbol*16
   TYPE(tpfun_expression), pointer :: exprot
   double precision, allocatable :: dummy(:)
   integer, allocatable :: idum(:)
   double precision tdummy
!---------------------
!>>>>> 20: no of TPFUNS
   read(lin)i,j,kx,freetpfun
!   write(*,*)'TPFUNS: ',i,j,kx,freetpfun
   if(i.ne.tpfun_expression_version) then
      write(*,*)'Wrong version of tpfun_expression'
      gx%bmperr=7777; goto 1000
   endif
   if(j.ne.tpfun_root_version) then
      write(*,*)'Wrong version of tpfun_root'
      gx%bmperr=7777; goto 1000
   endif
   if(kx.ne.tpfun_parres_version)then
      write(*,*)'Wrong version of tpfun_perres'
      gx%bmperr=7777; goto 1000
   endif
! tpfunction root list
   allfun: do j=1,freetpfun-1
! many symbols may already be entered, skip those without flagging error
!      read(lin)tpfuns(j)%symbol,tpfuns(j)%noofranges,tpfuns(j)%nextfree
!>>>>> 21: symbol and ranges
      read(lin)symbol,ranges,free
!      write(*,*)'Function: ',symbol,ranges,free
      call find_tpfun_by_name(symbol,lrot)
      if(gx%bmperr.eq.0 .and. skip) then
! already entered TPFUN, but we must read the rest anyway
!>>>>> 22A: dummy read number of ranges
!         write(*,*)'dummy read 1',j,ranges,freetpfun,symbol
         if(ranges.eq.0) then
! a "function" with just a value
            read(lin)tdummy
!            write(*,*)'dummy read 17: ',tdummy
            cycle allfun
         endif
         allocate(dummy(ranges))
         read(lin)(dummy(i),i=1,ranges)
         deallocate(dummy)
         do kx=1,ranges
!>>>>> 23A: dummy read number of coefficients in this range
!            write(*,*)'dummy read 2',j,kx
            read(lin)nc
            allocate(dummy(nc))
            allocate(idum(nc))
!>>>>> 24A: dummy read of coefficients, powers etc in this range
!            write(*,*)'dummy read 3',j,nc
            read(lin)(dummy(i),i=1,nc)
            read(lin)(idum(i),i=1,nc)
            read(lin)(idum(i),i=1,nc)
            read(lin)(idum(i),i=1,nc)
            read(lin)(idum(i),i=1,nc)
            read(lin)(idum(i),i=1,nc)
            deallocate(dummy)
            deallocate(idum)
         enddo
!>>>>> 25A: dummy read of high T limit
!         write(*,*)'dummy read 4',j
         read(lin)tdummy
!         write(*,*)'dummy read 5',j,tdummy
         cycle allfun
      elseif(gx%bmperr.eq.4060) then
! new TPFUN
!         write(*,*)'real read',j
         gx%bmperr=0;
         tpfuns(j)%symbol=symbol; tpfuns(j)%noofranges=ranges;
         tpfuns(j)%nextfree=free
      else
!         write(*,*)'error read',j,ranges,freetpfun,symbol
         goto 1000
      endif
! ----------------- read read
      nr=tpfuns(j)%noofranges
      if(nr.gt.0) then
         allocate(tpfuns(j)%limits(nr))
         allocate(tpfuns(j)%funlinks(nr))
!>>>>> 22B: read tlimits
         read(lin)(tpfuns(j)%limits(i),i=1,nr)
! read the expressions
         do kx=1,nr
            exprot=>tpfuns(j)%funlinks(kx)
!>>>>> 23B: dummy read number of coefficients in this range
            read(lin)exprot%noofcoeffs
            nc=exprot%noofcoeffs
            if(nc.gt.0) then
               allocate(exprot%coeffs(nc))
               allocate(exprot%tpow(nc))
               allocate(exprot%ppow(nc))
               allocate(exprot%wpow(nc))
               allocate(exprot%plevel(nc))
               allocate(exprot%link(nc))
!>>>>> 24A: read of coefficients, powers etc in this range
               read(lin)(exprot%coeffs(i),i=1,nc)
               read(lin)(exprot%tpow(i),i=1,nc)
               read(lin)(exprot%ppow(i),i=1,nc)
               read(lin)(exprot%wpow(i),i=1,nc)
               read(lin)(exprot%plevel(i),i=1,nc)
               read(lin)(exprot%link(i),i=1,nc)
            endif
         enddo
      endif
!>>>>> 25B: read of high T limit
      read(lin)tpfuns(j)%hightlimit
   enddo allfun
! We assume TPFUNS have a correct free list but we have to set freetpfun ???
!   freetpfun
!----------- tpfunction result list
!   read(lin)(tpres(j)%tpused,tpres(j)%results,j=1,freetpfun-1)
1000 continue
   return
 end subroutine tpfunread

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

! END MODULE TPFUNLIB

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>      20. Subroutines used by applications
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine makeoptvname(name,indx)
    implicit none
    character name*(*)
    integer indx
!\end{verbatim}
    if(indx.lt.99) then
       if(indx.le.9) then 
          name(1:2)='A0'
          name(3:3)=char(indx+ichar('0'))
       else
          name(1:1)='A'
          name(2:2)=char(indx/10+ichar('0'))
          name(3:3)=char(mod(indx,10)+ichar('0'))
       endif
    else
       name='A99'
    endif
1000 continue
    return
  end subroutine makeoptvname

