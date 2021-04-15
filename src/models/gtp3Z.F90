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
!
!\addtotable subroutine tpfun_init
!\begin{verbatim}
 subroutine tpfun_init(nf,tpres)
! allocate tpfuns and create a free list inside the tpfuns
   implicit none
   integer nf
! use tpres declared externally for parallel processing
  TYPE(tpfun_parres), dimension(:), allocatable :: tpres
!\end{verbatim}
   integer ifri
   allocate(tpfuns(nf))
! tpres allocated when creating equilibria
!   allocate(tpres(nf))
! create free list for named functions records
   freetpfun=1
   do ifri=1,nf-1
      tpfuns(ifri)%nextorsymbol=ifri+1
      tpfuns(ifri)%noofranges=0
      tpfuns(ifri)%status=0
      tpfuns(ifri)%forcenewcalc=0
! should also be initiallized ??
!      tpres(ifri)%forcenewcalc=0
   enddo
! The last TP function has no next link
   tpfuns(nf)%nextorsymbol=-1
   return
 END SUBROUTINE tpfun_init

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function notpf
!\begin{verbatim}
 integer function notpf()
! number of tpfunctions because freetpfun is private
   implicit none
!\end{verbatim}
   notpf=freetpfun-1
 end function notpf

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_tpfun_by_name
!\begin{verbatim}
 subroutine find_tpfun_by_name(name,lrot)
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

!\addtotable subroutine find_tpfun_by_name_exact
!\begin{verbatim} %-
 subroutine find_tpfun_by_name_exact(name,lrot,notent)
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

!\addtotable subroutine eval_tpfun
!\begin{verbatim}
 subroutine eval_tpfun(lrot,tpval,result,tpres)
!    subroutine eval_tpfun(lrot,tpval,symval,result)
! evaluate a TP function with several T ranges
   implicit none
   integer lrot
   double precision tpval(2),result(6),xxx
! changes to avoid memory leak in valgrind
   TYPE(tpfun_parres), dimension(*) :: tpres
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
   if(lrot.le.0) then
      result=zero
      goto 1000
   elseif(btest(tpfuns(lrot)%status,TPCONST)) then
! TP symbol is a constant, value stored in tpfuns(lrot)%limits(1)
! This takes care of updating assessment parameters!!
      result=zero
      result(1)=tpfuns(lrot)%limits(1)
! wow, we must not forget to store the constant in tpres(lrot)%results!
      goto 990
   else
! check if previous values can be used
! tpfuns(lrot)%forcenewcalc is located with the function expression
! tpres(lrot)%forcenewcalc is different for each ceq, there can be several
! IT IS MEANINGLESS TO COMPARE THEM ... 
      if(tpres(lrot)%forcenewcalc.eq.tpfuns(lrot)%forcenewcalc) then
         if(abs(tpres(lrot)%tpused(1)-tpval(1)).le.&
              mini*tpres(lrot)%tpused(1) .and. &
              (abs(tpres(lrot)%tpused(2)-tpval(2)).le.&
              mini*tpres(lrot)%tpused(2))) then
            result=tpres(lrot)%results
!            write(*,12)'3Z oldval: ',lrot,tpres(lrot)%forcenewcalc,&
!                 tpfuns(lrot)%forcenewcalc,tpres(lrot)%results(1),tpval(1)
!12          format(a,i5,2i4,4(1pe12.4))
            goto 1000
         endif
!      else
!         write(*,*)'3Z forced recalc: ',lrot,tpres(lrot)%forcenewcalc,&
!              tpfuns(lrot)%forcenewcalc
      endif
! new values must be calculated
!         write(*,23)'3Z new T,P: ',lrot,tpres(lrot)%tpused,tpval
!23       format(a,i4,4(1pe12.4))
!         result=zero
   endif
! we must calculate the function
!   write(*,35)'3Z new TPval:',lrot,tpfuns(lrot)%forcenewcalc,&
!        tpres(lrot)%forcenewcalc,&
!        abs(tpres(lrot)%tpused(1)-tpval(1)),abs(tpres(lrot)%tpused(2)-tpval(2))
!35 format(a,3i5,2(1pe12.4))
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
   if(gx%bmperr.ne.0) then
      write(*,901)gx%bmperr,tpfuns(lrot)%symbol
901   format('Error ',i5,' evaluating tp function: ',a)
      goto 1000
   endif
990 continue
!   new: do i=1,6
!      tpres(lrot)%results(i)=result(i)
!   enddo new
!   xxx=tpres(lrot)%results(1)
   tpres(lrot)%results=result
   tpres(lrot)%forcenewcalc=tpfuns(lrot)%forcenewcalc
   tpres(lrot)%tpused(1)=tpval(1)
   tpres(lrot)%tpused(2)=tpval(2)
! Searching for strange bug when entering parameter ...
!   write(*,991)'3Z new value: ',lrot,tpres(lrot)%forcenewcalc,&
!        tpres(lrot)%results(1),tpval(1),xxx
!991 format(a,2i5,6(1pe12.4))
!22  format(A,i3,4(1PE11.2))
1000 continue
   return
 end subroutine eval_tpfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine list_tpfun
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
! constant equal to zero ??
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
!20 format(a)
1000 continue
   return
 end subroutine list_tpfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine list_all_funs
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
!        ' Predefined symbols:'/&
!        ' BELOW(TB) = something;'/&
!        ' ABOVE(TB) = 1-BELOW(TB);'/&
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

!\addtotable subroutine list_unentered_funs
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
!1000 continue
   return
 end subroutine list_unentered_funs

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine ct1xfn
!\begin{verbatim}
 subroutine ct1xfn(string,ip,nc,coeff,koder,fromtdb)
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
! check consistency
   implicit none
   integer ip,nc,koder(5,*)
   character string*(*)
   double precision coeff(*)
   logical fromtdb
!\end{verbatim} %+
!   implicit double precision (a-h,o-z)
!   integer, parameter :: nunary=5
   integer, parameter :: nunary=6
   integer, parameter :: lenfnsym=16
   double precision, parameter :: zero=0.0D0,one=1.0D0
   integer i,j,jss,levelp,mterm,ipower,nterm
   double precision sign,val,another
   character ch1*1
   logical zeroc
   character symbol*(lenfnsym),unary(nunary)*6
   character*2, parameter :: tsym='T ',psym='P '
! NEIN is the Einstein function
! MAX1 is 1.0 if argument is larger than 1.0, error if argument negative
! LOG is LOG10 and LN is the natural logarithm!!!
   DATA unary/'LOG   ','LN    ','EXP   ','ERF   ','GEIN','MAX1  '/
!   DATA unary/'LOG   ','LN    ','EXP   ','ERF   ','INTEIN','MAX1  '/
!   DATA unary/'LOG   ','LN    ','EXP   ','ERF   ','XNEIN ','MAX1  '/
!   DATA unary/'LOG   ','LN    ','EXP   ','ERF   ','XNEIN '/
!   DATA unary/'LOG   ','LN    ','EXP   ','ERF   ','ABOVE ','BELOW '/
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
      write(*,*)'ct1xfn 66:',ip,ch1,' >',trim(string),'<'
      gx%bmperr=4001
      goto 1000
   endif
!   write(*,202)ip,string(1:ip+5)
!202 format('Expected real at position: ',i5,' in >',a,'< ')
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
      call store_tpfun_dummy(symbol)
   else
! otherwise give error message
      write(*,*)'TPFUN contain unknown symbol: ',symbol,freetpfun-1
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
!      write(*,*)'3Z we found a multiplication: ',trim(string),ip
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
      elseif(ch1.ge.'0' .and. ch1.le.'9') then
! multiplying value in coeff(nterm)  with another number
!         write(*,*)'3Z string and position: "',trim(string),'"',ip
!         write(*,*)'3Z multplication followed by digit: ',ch1,nterm
!         write(*,*)'3Z data: ',buperr,coeff(nterm),trim(string(ip:))
! Does getrel increment ip?? NO
         call getrel(string,ip,another)
         if(buperr.ne.0) then
!            write(*,*)'2Z error from getrel',buperr
            gx%bmperr=buperr; goto 1000
         endif
!         write(*,'(a,i2,2(1pe12.4))')'3Z multiplying two numbers: ',nterm,&
!              coeff(nterm),another
         coeff(nterm)=coeff(nterm)*another
! now we expect an operator or ) or ;
         goto 400
!      else
! multiply symbol with something ....
!         write(*,*)'ct1xfn 4: ',nterm,ip,coeff(nterm)
!654      format(a,i5,'"',a,'"',i5)
      endif
! new we expect a symbol or end of expression
      goto 200
   elseif(ch1.eq.';') then
      goto 900
   endif after
   write(*,777)ch1,ip,trim(string)
777 format('3Z Illegal character "',a,'" at pos ',i3,' in expression "',a,'"')
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
!990 format('ct1xfn 99> ',1PE15.6,5I7)
1000 continue
   return
 end subroutine ct1xfn

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine ct1getsym
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
!1000 return
 end subroutine ct1getsym

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine ct1power
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

!\addtotable subroutine ct1mfn
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
! special for unformatted files, lrot < 0 and this index MUST be used
! ignore freetpfun!!
   integer ir
   character name*16
   if(lrot.lt.0) then
! store funtion at this specific place!!
      lrot=-lrot
!      if(lrot.gt.freetpfun) then
         write(*,*)'Storing at position above freetpfun',lrot
!         gx%bmperr=4399; goto 1000
!      endif
   else
      lrot=freetpfun
!   write(*,*)'ct1mfn: ',freetpfun
!   write(*,*)'ct1mfn: ',lrot,tpfuns(lrot)%nextorsymbol
      if(lrot.gt.0) then
         freetpfun=tpfuns(lrot)%nextorsymbol
         tpfuns(lrot)%nextorsymbol=0
      else
! no more tpfun records
         write(*,*)'No more space for TP functions: ',size(tpfuns)
         gx%bmperr=4014; goto 1000
      endif
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

!\addtotable subroutine ct2mfn
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
         gx%bmperr=4348; goto 1000
      endif
   else
! illegal value of lrot
      gx%bmperr=4349; goto 1000
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

!\addtotable subroutine ct1mexpr
!\begin{verbatim} %-
 subroutine ct1mexpr(nc,coeff,koder,lrot)
!...makes a datastructure of an expression. root is returned in lrot
!   implicit double precision (a-h,o-z)
   implicit none
   integer nc,koder(5,*)
!   TYPE(tpfun_expression), pointer :: lrot
   TYPE(tpfun_expression) :: lrot
!   TYPE(tpfun_expression), pointer :: noexpr
   double precision coeff(*)
!\end{verbatim} %+
   integer i
!   write(*,*)'3Z in ct1mexpr',nc
   lrot%noofcoeffs=nc
   if(nc.le.0) then
!      nullify(lrot)
      goto 1000
   endif
! allocate an expression record and then allocate all arrays
!   allocate(lrot)
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
1000  continue
!   write(*,*)'3Z leaving ct1mexpr'
   return
 end subroutine ct1mexpr

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine ct1efn
!\begin{verbatim} %-
 subroutine ct1efn(inrot,tpval,val,tpres)
!...evaluates a datastructure of an expression. Value returned in val
!     inrot is root expression tpfunction record
!     tpval is valuse of T and P,
!     val is array of values calculated here
!     tpres is array of all calculated functions
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
   TYPE(tpfun_parres), dimension(*) :: tpres
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
!   write(*,*)'ct1efn 0: ',tpval(1)
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
!...power of symbols is handeled below
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
! nonzero link4 inserted in the wrong term in step3.OCM ...
      link4=0
!       if(link.ne.0) write(*,201)'funev: ',lrot,ic,link,link3
!201    format(a,10i5)
!      write(*,'(a,5i4,1pe12.4)')'3Z intein2: ',ic,ipow,link,link3,link4,ff
      if(link.lt.0 .and. link3.gt.1000) then
! if link is negative (unary funktion) and link3 is >1000 (link)
! we must evaluate link3 first
         link4=link3-1000
         if(abs(tpres(link4)%tpused(1)-tpval(1)).lt.&
              mini*tpres(link4)%tpused(1) .and. &
              abs(tpres(link4)%tpused(2)-tpval(2)).lt.&
              mini*tpres(link4)%tpused(2)) then
! The test for forcenewcalc is not reasonable:
! tpres%forcenewcalc is local for each equilibrium
! tpfun(link4)%forcenewcalc is global for whole system.  If we calculate in
!   parallel there is no reason they should be the same
! It creates problem in testcond1.OCM after adding a mobility parameter
!   although I do not understand why that should create the problem
! I added it to be sure that updated assessment parameters should be used
! but there seems no problem with that ...
! Better to speed a few months to rerwrite the whole TPFUN package ...
!              mini*tpres(link4)%tpused(2) .and. &
! added test of forcenewcalc ... removed ??
!              tpres(link4)%forcenewcalc.eq.tpfuns(link4)%forcenewcalc) then
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
            tpres(link4)%forcenewcalc=tpfuns(link4)%forcenewcalc
! DANGER changing wpow
!            write(*,22)'3Z wpow 1:      ',level,nc,ic,exprot%wpow(ic),link-1000
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
! DANGER, restoring wpow, note also wpow pop below!
!         write(*,22)'3Z wpow 2:      ',level,nc,ic,exprot%wpow(ic),link4
!22       format(a,7i7)
         exprot%wpow(ic)=link4
      endif
!-------------------------------------------------------------
!      write(*,'(a,5i4,1pe12.4)')'3Z intein3: ',ic,ipow,link,link3,link4,ff
      evlink: if(link.gt.0) then
! link to another symbol, extract its value and use chain rule
! extract the results from the symbol if already calculated
! if not calculated then do that and then recalculate this term
         linkif: if(abs(tpres(link)%tpused(1)-tpval(1)).lt.&
              mini*tpres(link)%tpused(1) .and. &
              abs(tpres(link)%tpused(2)-tpval(2)).lt.&
!              mini*tpres(link)%tpused(2)) then
              mini*tpres(link)%tpused(2) .and. &
! added this check as it seems new assessment coefficients are nor used!!
              tpres(link)%forcenewcalc.eq.tpfuns(link)%forcenewcalc) then
! Valgrid complained about uninitial variable in if above, I do not know which
            jpow=exprot%wpow(ic)
!---------------------------------------------
            jpowif: if(jpow.gt.1000) then
! suck, two functions have to be multiplied ....
               link2=jpow-1000
               jpowev: if(abs(tpres(link2)%tpused(1)-tpval(1)).lt.&
                    mini*tpres(link2)%tpused(1) .and. &
                    abs(tpres(link2)%tpused(2)-tpval(2)).lt.&
!                    mini*tpres(link2)%tpused(2)) then
                    mini*tpres(link2)%tpused(2) .and. &
! added this check as it seems new assessment coefficients are nor used!!
                  tpres(link2)%forcenewcalc.eq.tpfuns(link2)%forcenewcalc) then
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
                  tpres(link2)%forcenewcalc=tpfuns(link2)%forcenewcalc
                  if(btest(tpfuns(link2)%status,TPCONST)) then
!                     write(*,*)'3Z Link to a constant 1',&
!                          link2,tpfuns(link2)%limits(1)
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
! MEMORY LEAK
                  nullify(temp)
                  topsave%exprot=>exprot
                  topsave%level=level
                  topsave%saveic=ic-1; topsave%savenc=nc
!                  write(*,22)'3Z wpow push 1: ',level,nc,ic,0,link4
                  topsave%savetp=link2; topsave%savelink4=link4
                  topsave%saveval=val
                  link4=0
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
            tpres(link)%forcenewcalc=tpfuns(link)%forcenewcalc
            if(btest(tpfuns(link)%status,TPCONST)) then
! the function is a constant!!
!               write(*,*)'3Z Link to a constant 2',link,tpfuns(link)%limits(1)
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
! MEMORY LEAK
            nullify(temp)
            topsave%exprot=>exprot
            topsave%level=level
            topsave%saveic=ic-1; topsave%savenc=nc
!            write(*,22)'3Z wpow push 2: ',level,nc,ic,0,link4
            topsave%savetp=link; topsave%savelink4=link4
            topsave%saveval=val
            link4=0
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
!         write(*,'(a,5i4,1pe12.4)')'3Z intein4: ',ic,ipow,link,link3,link4,ff
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
            write(*,*)'TP ipow error: ',ipow
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
                 mini*tpres(link2)%tpused(2) .and. &
! added this check as it seems new assessment coefficients are nor used!!
                 tpres(link2)%forcenewcalc.eq.tpfuns(link2)%forcenewcalc) then
               symval=tpres(link2)%results
            else
! one must evaluaste another function, it is recursive through eval_tpfun
! which will call ct1efn again but this is handelled automatically?????
! One should add some check that two TP functions does not call each other
! to infinite depth
               tpres(link2)%forcenewcalc=tpfuns(link2)%forcenewcalc
               if(btest(tpfuns(link2)%status,TPCONST)) then
! the function is a constant!!
!                  write(*,*)'3Z link to a constant 3',link2,&
!                       tpfuns(link2)%limits(1)
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
! MEMORY LEAK
               nullify(temp)
               topsave%exprot=>exprot
               topsave%level=level
               topsave%saveic=ic-1; topsave%savenc=nc
! this searching for strange bug at midsummer 2018 ...
!               write(*,22)'3Z wpow push 3: ',level,nc,ic,0,link4
               topsave%savetp=link2; topsave%savelink4=link4
               topsave%saveval=val
               link4=0
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
! -5: GEIN, is the Einstein function, integrated as a Gibbs energy
!             the argument is the Einstein T
! -6: MAX1, if argument <0 ERROR, if >1 replace by 1
!         write(*,'(a,5i4,1pe12.4)')'3Z intein5: ',ic,ipow,link,link3,link4,ff
         evunfun: if(unfun.eq.-1) then
! LOG base 10
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
! LN NATURAL LOGARITHM
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
! EXPonential
! ff=ff*exp(gg)
            d2fdp2=exp(gg)*(d2fdp2+2.0D0*dfdp*dgdp+ff*d2gdp2+ff*(dgdp)**2)
            d2fdtdp=exp(gg)*(d2fdtdp+dfdt*dgdp+dfdp*dgdt+ff*d2gdtdp &
                 +ff*dgdt*dgdp)
            d2fdt2=exp(gg)*(d2fdt2+2.0D0*dfdt*dgdt+ff*d2gdt2+ff*(dgdt)**2)
            dfdp=exp(gg)*(dfdp+ff*dgdp)
            dfdt=exp(gg)*(dfdt+ff*dgdt)
            ff=ff*exp(gg)
         elseif(unfun.eq.-4) then
! ERROR FUNCTION or ABOVE not implemented
            write(*,*)'Error function not implemented'
            stop 71
         elseif(unfun.eq.-5) then
! INTEGRATED EINSTEIN: GEIN = 1.5*R*THETA + 3*R*T*LN(EXP(THETA/T)+1), THETA=gg
            if(dfdt.ne.zero) then
               write(*,*)'3Z GEIN must not be multiplied with T!'
               gx%bmperr=4399; goto 1000
            endif
!           write(*,'(a,5i5,1pe12.4)')'3Z intein6: ',ic,ipow,link,link3,link4,ff
! ff is the constant factor in front of the Einstein function
! It is overwritten by the Einsten function (multiplied by original ff)
            call tpfun_geinstein(tpval,gg,ff,dfdt,dfdp,d2fdt2,d2fdtdp,d2fdp2)
!            write(*,'(a,i3,6(1pe12.4))')'3Z call Einstein:',link,&
!                 gg,ff,dfdt,d2fdt2
! ff is the coefficient for the Einstein Functions, should be a constant ...?
!            write(*,*)'Einstein function not implemented'
!            stop 72
            if(gx%bmperr.ne.0) goto 1000
         elseif(unfun.eq.-6) then
! MAX1 function, used for SRO .... function and derivatives in gg, dgdt etc.
!            write(*,*)'MAX1 function',gg
            if(gg.le.zero) then
               write(*,*)'MAX1 called with negative argument',gg
               stop 73
            endif
            if(gg.le.one) then
! just copy values from g to f
               d2fdp2=d2gdp2; d2gdtdp=d2fdtdp; d2fdt2=d2gdt2
               dfdp=dgdp; dfdt=dgdt; ff=gg
            else
! function value is 1 and all derivatives zero
               d2fdp2=zero; d2fdtdp=zero; d2fdt2=zero
               dfdp=zero; dfdt=zero; ff=one
            endif
         else
! undefined function
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
! For some unknown reason topsave%saveic is ic-1 !!! correct below
!      write(*,22)'3Z wpow pop 1:  ',level,nc,ic,0,link4
      exprot=>topsave%exprot
! MEMORY LEAK avoided by deallocate topsave ??
!      write(*,*)'Trying to remove memory leak'
      temp=>topsave%previous
      deallocate(topsave)
!      write(*,*)'Deallocated topsave'
!      topsave=>topsave%previous
      topsave=>temp
      level=level-1
! restart from coefficient ic, note the value saved is ic-1 !!
      if(ic.ge.0 .and. ic.lt.nc) then
! restore value in %wpow !!!
! without this an expression like VCRBCC*EXP(ZCRBCC) became
! just EXP(ZRBCC) as the link to VCRBCC had been removed ...
! BUT for macro step3 the link4 was inserted in the wrong term !!!
!         if(link4.gt.1000 .and. exprot%wpow(ic).lt.1000) then
         if(link4.gt.1000 .and. exprot%wpow(ic+1).lt.1000) then
!            write(*,22)'3Z wpow save:   ',level+1,nc,ic,exprot%wpow(ic),link4
! topsave%saveic is ic-1, I do not know why but this correction is added now!
            exprot%wpow(ic+1)=link4
         endif
      endif
      goto 200
    endif
!
1000 continue
   return
 end subroutine ct1efn !level %wpow link4

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine tpfun_geinstein
!\begin{verbatim}
 subroutine tpfun_geinstein(tpval,gg,ff,dfdt,dfdp,d2fdt2,d2fdtdp,d2fdp2)
! evaluates the integrated Einstein function (including 1.5*R) INTEIN/GEIN
! gg is the value of the Einstein THETA
! ff is a constant factor which should be multiplied with all terms
! ff is overwritten with the Einstein function (multiplied with ff in)
! the other parameters are derivatives of the integrated Einstein function
   implicit none
   double precision tpval(*)
   double precision gg,ff,dfdt,dfdp,d2fdt2,d2fdtdp,d2fdp2
!\end{verbatim}
   double precision kvot,kvotexpkvotm1,expmkvot,lnexpkvot,ww,rgas
! return ff = 1.5*R*gg + 3*R*T*LN(EXP(-gg/T) + 1) and derivatives
! gg must be a constant >0
   rgas=globaldata%rgas
!   write(*,*)'3Z in Einstein function',gg,tpval(1)
   ww=ff
   kvot=gg/tpval(1)
   if(kvot.gt.2.0d2) then
! handle extreme values of kvot, we divide by kvotexpkvotm1**2 by expmkvot bolw
      expmkvot=one
      kvotexpkvotm1=zero
      lnexpkvot=zero
!      write(*,'(a,5(1pe12.4))')'3Z Einetsin 1: ',kvot,expmkvot,&
!           kvotexpkvotm1,lnexpkvot
   else
      expmkvot=exp(-kvot)
      kvotexpkvotm1=kvot/(exp(kvot)-one)
      lnexpkvot=log(one-expmkvot)
!      write(*,'(a,5(1pe12.4))')'3Z Einetsin 2: ',kvot,expmkvot,&
!           kvotexpkvotm1,lnexpkvot
   endif
! this is the integral G contribution from an Einstein solid
   ff=1.5d0*rgas*gg*ww + 3.0D0*rgas*tpval(1)*lnexpkvot*ww
   dfdt=3.0d0*rgas*(lnexpkvot-kvotexpkvotm1)*ww
!   write(*,10)rgas,kvot,lnexpkvot,kvotexpkvotm1,dfdt
!10 format('3Z bug: ',6(1pe12.4))
   dfdp=zero
! this is the second derivative of G wrt T; i.e. the Einstein solid Cp equation
   d2fdt2=-3.0d0*rgas*kvotexpkvotm1**2/(expmkvot*tpval(1))*ww
   d2fdtdp=zero
   d2fdp2=zero
1000 continue
   return
 end subroutine tpfun_geinstein

   !/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine ct1wfn
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
! these should be the same as in ct1xfn !!! ??
   DATA unary/'LOG   ','LN    ','EXP   ','ERF   ','GEIN  ','MAX1  '/
!   DATA unary/'LOG   ','LN    ','EXP   ','ERF   ','INTEIN','MAX1  '/
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
!71    format(A,I5,1PE15.6,5I5)
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
!230    continue
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
! if 4th argument >0 then write a sign
      call wrinum(string,ip,12,6,coeff(level))
      mult=-1
   else
! in the case of a single value exactly 1 without unary or T or P power
! the number was never written
!      write(*,203)'ct1wfn2: ',(koder(i,level),i=1,4),coeff(level)
!203   format(a,4i4,1pe12.4)
      do i=1,4
         if(koder(i,level).ne.0) goto 219
      enddo
! without this the Inden magnetic function will miss its initial 1.0
!      call wrinum(string,ip,2,0,coeff(level))
! changed 20.03.17/BoS because EXP(T)+1 missed the + between ) and 1
! Force wrinum to write positive signs by 4th parameter positive
      call wrinum(string,ip,2,1,coeff(level))
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

!\addtotable subroutine ct1wpow
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

!\addtotable subroutine enter_tpfun_interactivly
!\begin{verbatim}
 subroutine enter_tpfun_interactivly(cline,ip,longline,jp)
! interactive input of a TP expression, whole function returned in longline
!   implicit double precision (a-h,o-z)
   implicit none
   integer ip,jp
   character cline*(*),longline*(*)
!\end{verbatim}
   character line*80,ch1*1
   integer nexpr,lsc,kkp
   double precision xx
!   write(*,*)'Max ',len(longline),' characters'
   call gparrdx('Low temperature limit: ',cline,ip,xx,2.9815D2,'?Enter TPfun')
   if(buperr.ne.0) then
! set default low limit
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
   call gparcx('Give expression, end with ";":',cline,ip,6,line,';',&
        '?Enter TPfun')
   if(buperr.ne.0) then
      buperr=0; line=';'
   endif
120 continue
   longline(jp:)=line
   jp=len_trim(longline)+1
!   write(*,*)'tpfun: ',longline(1:jp)
! lsc is position after the ";" in any previous range
   if(index(longline(lsc:),';').le.0) then
      call gparcx('&',cline,ip,6,line,';','?Enter TPfun')
      if(buperr.ne.0) then
         buperr=0; line=';'
      endif
      goto 120
   endif
!150 continue
! make sure there is a ; at the end of each expression
   kkp=index(longline(nexpr:),';')
!   write(*,130)'3Z pos1: ',nexpr,kkp,lsc,jp,trim(longline)
!130 format(a,4i4,': ',a/26x,'123456789.123456789.123456789.123456789.')
!   write(*,*)'tpfun add ;'
   if(kkp.le.0) then
      kkp=len_trim(longline)
      longline(kkp+1:)='; '
      jp=kkp+3
      nexpr=jp
      write(*,*)'3Z adding ; at position: ',kkp+1,nexpr
   else
!      nexpr=kkp+1
      nexpr=len_trim(longline)+2
   endif
! lsc is position of ; for previous range
!   write(*,130)'3Z pos2: ',nexpr,kkp,lsc,jp,trim(longline)
   lsc=nexpr
   call gparrdx('Upper temperature limit ',cline,ip,xx,6.0D3,'?Enter TPfun')
   if(buperr.ne.0) then
      buperr=0; xx=6.0D3
   endif
! enter a space after ;
   jp=jp+1
   call wrinum(longline,jp,8,0,xx)
   if(buperr.ne.0) goto 1000
   call gparcdx('Any more ranges',cline,ip,1,ch1,'N','?Enter TPfun')
!   write(*,*)'3Z ch1: ',ch1
   if(ch1.eq.'n' .or. ch1.eq.'N') then
      longline(jp:)=' N'
      jp=jp+3
   else
      longline(jp:)=' Y'
      jp=jp+3
      goto 115
   endif
! remove any "#" (comes from TC functions)
900 continue
   kkp=index(longline,'#')
   if(kkp.gt.0) then
      longline(kkp:kkp)=' '
      goto 900
   endif
!   write(*,910)'3Z tpf: ',jp,trim(longline)
!910 format(a,i3,': ',a)
!
1000 continue
   return
 end subroutine enter_tpfun_interactivly

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine tpfun_deallocate
!\begin{verbatim}
 subroutine tpfun_deallocate
! deallocates all arrays associated with a TP function
!\end{verbatim}
   implicit none
   TYPE(tpfun_expression), pointer :: exprot
   integer j,nr,nexp,nc
!   write(*,*)'3Z freetpfun: ',freetpfun
   do j=1,freetpfun-1
      nr=tpfuns(j)%noofranges
      if(nr.gt.0) then
! modified 170517 due to memory leaks when read/write unformatted
         do nc=1,nr
            exprot=>tpfuns(j)%funlinks(nc)
!            write(*,*)'3Z Deallocating TP function',j,nc
            deallocate(exprot%tpow)
            deallocate(exprot%ppow)
            deallocate(exprot%wpow)
            deallocate(exprot%plevel)
            deallocate(exprot%link)
            deallocate(exprot%coeffs)
         enddo
!
         deallocate(tpfuns(j)%funlinks)
         deallocate(tpfuns(j)%limits)
      endif
   enddo
   deallocate(tpfuns)
!1000 continue
   return
 end subroutine tpfun_deallocate

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine store_tpfun_dummy
!\begin{verbatim}
 subroutine store_tpfun_dummy(symbol)
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
      freetpfun=tpfuns(lrot)%nextorsymbol
      tpfuns(lrot)%nextorsymbol=0
   else
      write(*,*)'No space for TP functions: ',size(tpfuns)
      gx%bmperr=4014; goto 1000
   endif
   tpfuns(lrot)%noofranges=0
   name=symbol
   call capson(name)
   tpfuns(lrot)%symbol=name
   tpfuns(lrot)%status=ibset(tpfuns(lrot)%status,TPNOTENT)
   tpfuns(lrot)%rewind=0
1000 continue
   return
 end subroutine store_tpfun_dummy

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine store_tpfun
!\begin{verbatim}
 subroutine store_tpfun(symbol,text,lrot,rewind)
! creates a data structure for a TP function called symbol with several ranges
! text is whole expression
! lrot is returned as index.  If fromtdb is FALSE and lrot<0 it is a new
!                             expression for an old symbol
! if fromtdb is TRUE references to unknown functions are allowed
! default low temperature limit is 298.16; high 6000
   implicit none
   integer lrot,rewind
   character*(*) text,symbol
!   logical fromtdb
!\end{verbatim}
! max number of ranges, max number of coefficents in each range
!   integer, parameter :: mrange=20,mc=15
! in a paper more than 15 terms were used for a TP function!
   integer, parameter :: mrange=20,mc=20
   integer jss,nc,ip,nrange,cbug
   real tlim(mrange)
   double precision coeff(mc),val
   integer koder(5,mc)
! attempt to remove big memory leak
!   TYPE(tpfun_expression) :: links(mrange)
!   TYPE(tpfun_expression), target :: links(mrange)
   TYPE(tpfun_expression) :: links(mrange)
!   TYPE(tpfun_expression), pointer :: ltpexpr
   character ch1*1,lsym*(lenfnsym)
   logical already,fromtdb
! check if function already entered, there are freetpfun-1 of them
! ignore functions that start with a "_" as they are parameters
!   lrot=0
! special when read unformatted or direct files, lrot<0 and this ??
! must be the location for storing the function ...
   fromtdb=.TRUE.
   if(rewind.lt.0) fromtdb=.FALSE.
   already=.FALSE.
   if(symbol(1:1).ne.'_') then
      lsym=symbol
      call capson(lsym)
!      write(*,*)'3Z store_tpfun: ',trim(lsym),lrot,rewind
      do jss=1,freetpfun-1
!         write(*,17)jss,lsym,tpfuns(jss)%symbol
!17       format('enter_tpfun: ',i5,' >,',a,'=',a,'?')
         if(lsym.eq.tpfuns(jss)%symbol) then
            if(btest(tpfuns(jss)%status,TPNOTENT)) then
! function name already entered, now enter expression, this is from TDB files
               lrot=jss; already=.TRUE.
! mark the expression was entered at current rewind
               tpfuns(jss)%rewind=rewind; goto 18
            else
!               write(*,*)'amend tpfun? ',trim(lsym),fromtdb,lrot
               if(.NOT.fromtdb .and. lrot.lt.0) then
! this is an AMEND TPFUN, delete old expression to be able to store a new
                  lrot=jss; already=.TRUE.
                  nrange=tpfuns(lrot)%noofranges
!                  write(*,*)'Deallocating: ',lrot,nrange
                  deallocate(tpfuns(lrot)%limits)
                  deallocate(tpfuns(lrot)%funlinks)
                  tpfuns(lrot)%noofranges=0
                  tpfuns(lrot)%status=ibset(tpfuns(lrot)%status,TPNOTENT)
! we should clear the stored values! But those are stored separatly in all ceq
                  nrange=0; goto 18
               else
                  write(*,*)'3Z A never never error again! ',trim(symbol)
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
   cbug=ip
   call getrel(text,ip,val)
   if(buperr.ne.0) then
! A , has been used to select default low temperature limit
      if(text(ip:ip).eq.',') then
         buperr=0; val=298.15D0
      else
         write(*,*)'Illegal character for low temperature limit: ',text(ip:ip)
         val=298.15; buperr=0
!         write(*,19)ip,cbug,trim(text)
!19    format('TPFUN: ',2i3,' >',a)
      endif
! increement ip!
      ip=ip+1
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
! big memory leak ... still there ...
!      call ct1mexpr(nc,coeff,koder,ltpexpr)
!      links(nrange)=ltpexpr
!      ltpexpr=>links(nrange)
!      call ct1mexpr(nc,coeff,koder,ltpexpr)
!      write(*,*)'3Z calling ct1mexpr', nrange
      call ct1mexpr(nc,coeff,koder,links(nrange))
! attempt to remove memory leak
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
! force functions to be recalculated
!   write(*,*)'3Z calling force_recalculate from enter_tpfun'
   call force_recalculate_tpfuns
1000 continue
   return
 end subroutine store_tpfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine nested_tpfun
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
         gx%bmperr=4350; goto 1000
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

!\addtotable logical function compare_abbrev
!\begin{verbatim} %-
 logical function compare_abbrev(name1,name2)
! returns TRUE if name1 is an abbreviation of name2
! termintaes when a space is found in name1
! each part between _ or - can be abbreviated from the left
! a slash is treated as _
! case insensitive. Only 36 first characters compared
   implicit none
   character*(*) name1,name2   
!\end{verbatim}
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
! in species - has a special meaning of charge (and other things)
      if(ch1.eq.'-') then
         if(ch1.eq.lname2(jp:jp)) goto 300
         ch1='_'
      endif
! in species / has a special meaning of cluster
!      if(ch1.eq.'/') then
!         if(ch1.eq.lname2(jp:jp)) goto 300
!         write(*,*)'3Z accepting /'
!      endif
      if(ch1.eq.lname2(jp:jp)) goto 300
!      if(ch1.eq.'_' .or. ch1.eq.'-') then
      if(ch1.eq.'_') then
! we can abbreviate up to "_" in full name
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
!310    continue
   enddo bigloop
900 continue
   compare_abbrev=.TRUE.
1000 continue
   return
 end function compare_abbrev

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine enter_optvars
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
         freetpfun=tpfuns(lrot)%nextorsymbol
         tpfuns(lrot)%nextorsymbol=0
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

!\addtotable subroutine find_tpsymbol
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
   gx%bmperr=4351
   type=-1
200 continue
!1000 continue
   return
 end subroutine find_tpsymbol

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine store_tpconstant
!\begin{verbatim} %-
 subroutine store_tpconstant(symbol,value)
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
      freetpfun=tpfuns(lrot)%nextorsymbol
      tpfuns(lrot)%nextorsymbol=0
   endif
   allocate(tpfuns(lrot)%limits(1))
   allocate(tpfuns(lrot)%funlinks(1))
   call capson(symbol)
   tpfuns(lrot)%symbol=symbol
! mark this is a single value
   tpfuns(lrot)%status=ibset(tpfuns(lrot)%status,TPCONST)
200 continue
   tpfuns(lrot)%limits(1)=value
   nullify(tpfuns(lrot)%funlinks)
1000 continue
   return
 end subroutine store_tpconstant

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine change_optcoeff
!\begin{verbatim} %-
 subroutine change_optcoeff(lrot,value)
! change value of optimizing coefficient.  lrot is index
! -1 means just force recalculate
   implicit none
   integer lrot
   double precision value
!\end{verbatim} %+
   integer mrot
   if(lrot.gt.0 .and. lrot.lt.freetpfun-1) then
      if(.not.btest(tpfuns(lrot)%status,TPOPTCON)) then
         write(*,*)'Attempt to change non-existing coefficent',lrot
         gx%bmperr=7777; goto 1000
      endif
      tpfuns(lrot)%limits(1)=value
   endif
! force recalculation of all functions. HOW? the force_... does not work ...
   call force_recalculate_tpfuns
1000 continue
   return
 end subroutine change_optcoeff

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine force_recalculate_tpfuns
!\begin{verbatim} %-
 subroutine force_recalculate_tpfuns
! force recalculation of all tpfuns by incrementing an integer in tpfuns
!\end{verbatim} %+
   implicit none
   integer mrot
! it seems difficult to force recalculating all TP functions !!!
!   write(*,*)'3Z GLAVESCUMG: ',tpfuns(125)%forcenewcalc
   do mrot=1,freetpfun-1
      tpfuns(mrot)%forcenewcalc=tpfuns(mrot)%forcenewcalc+1
! I have no access to tpres here so I cannot see any current value ...
   enddo
!   write(*,*)'3Z Force recalculate tpfuns: ',freetpfun-1
!   write(*,*)'3Z GLAVESCUMG: ',tpfuns(125)%forcenewcalc
   return
 end subroutine force_recalculate_tpfuns

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_value_of_constant_name
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
!1000 continue
   return
 end subroutine get_value_of_constant_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_value_of_constant_index
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
!1000 continue
   return
 end subroutine get_value_of_constant_index

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_all_opt_coeff
!\begin{verbatim} %-
 subroutine get_all_opt_coeff(values)
! get values of all optimizing coefficients
   implicit none
   double precision values(*)
!\end{verbatim} %+
   write(*,*)'Not yet implemeneted'
!1000 continue
   return
 end subroutine get_all_opt_coeff

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine delete_all_tpfuns
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

!\addtotable subroutine save0tpfun
!\begin{verbatim}
 subroutine save0tpfun(lfun,iws,jfun)
! save one tpfun (or parameter) with index jfun in workspace iws
!   implicit double precision (a-h,o-z)
   implicit none
   integer lfun,iws(*),jfun
!\end{verbatim} %+
   integer nr,i,kx,nc,displace,lexpr,rsize,mmz
   TYPE(tpfun_expression), pointer :: exprot
   double precision dummy,xxx
! jfun can be zero meaning a parameter that is zero   
! unformatted
   if(jfun.eq.0) then
      iws(lfun)=0
   else
      nr=tpfuns(jfun)%noofranges
      rsize=4+nwch(16)+nr*(1+nwpr)+nwpr
      call wtake(lfun,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'Error reserving record for TPfun',buperr,rsize,nr
         gx%bmperr=4399; goto 1000
      endif
!      write(*,11)'3Z tpfun ',lfun,jfun,trim(tpfuns(jfun)%symbol),&
!           tpfuns(jfun)%noofranges,tpfuns(jfun)%status
!11    format(a,2i7,2x,a,2x,5i7)
      iws(lfun+1)=tpfuns(jfun)%noofranges
      iws(lfun+2)=tpfuns(jfun)%status
! what is nextfree??
      iws(lfun+3)=tpfuns(jfun)%nextorsymbol
      call storc(lfun+4,iws,tpfuns(jfun)%symbol)
      displace=4+nwch(16)
      call storrn(nr,iws(lfun+displace),tpfuns(jfun)%limits)
!         write(lut)(tpfuns(jfun)%limits(i),i=1,nr)
      call storr(lfun+displace+nr*nwpr,iws,tpfuns(jfun)%hightlimit)
! store location of expressions from displace
      displace=displace+nwpr*(nr+1)
! now the expressions, number of coefficients, nc, can be different
! link them from lfun
      do kx=1,nr
         exprot=>tpfuns(jfun)%funlinks(kx)
         nc=exprot%noofcoeffs
         rsize=1+nc*(5+nwpr)
         call wtake(lexpr,rsize,iws)
         if(buperr.ne.0) then
            write(*,*)'Error reserving record for TPfun'
            gx%bmperr=4399; goto 1000
         endif
         iws(lfun+displace+kx-1)=lexpr
         iws(lexpr)=nc
         mmz=lexpr+1
! The coefficients and the codes
         do i=1,nc
            iws(mmz)=exprot%link(i)
            iws(mmz+1)=exprot%tpow(i)
            iws(mmz+2)=exprot%ppow(i)
            iws(mmz+3)=exprot%wpow(i)
            iws(mmz+4)=exprot%plevel(i)
            call storr(mmz+5,iws,exprot%coeffs(i))
            call loadr(mmz+5,iws,xxx)
            mmz=mmz+5+nwpr
         enddo
      enddo
   endif
1000 continue
   return
 end subroutine save0tpfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine read0tpfun
!\begin{verbatim} %-
 subroutine read0tpfun(lfun,iws,jfun)
! read one TPfun from workspace
   implicit none
   integer lfun,jfun,iws(*)
!\end{verbatim}
   integer i,i2,kx,nc,nr,displace,lexpr,loklexpr,mmz
   TYPE(tpfun_expression), pointer :: exprot
   character*16 symbol
   double precision dummy
! jfun can be sero meaning no link to a TPFUN
!   read(lut)jfun,symbol,nr,i2
! the TPfuns are stored in an array, no need to allocate
!   if(iws(lfun).gt.0) then
   if(jfun.gt.0) then
      nr=iws(lfun+1)
      tpfuns(jfun)%noofranges=nr
      tpfuns(jfun)%status=iws(lfun+2)
      tpfuns(jfun)%nextorsymbol=iws(lfun+3)
      call loadc(lfun+4,iws,tpfuns(jfun)%symbol)
   else
      write(*,*)'not a function: ',lfun,jfun
      goto 1000
   endif
! special for optimizing variables
   if(btest(tpfuns(jfun)%status,TPOPTCON)) then
!      write(*,*)'3Z allocating zero limit for ',tpfuns(jfun)%symbol
      allocate(tpfuns(jfun)%limits(1))
      tpfuns(jfun)%limits(1)=zero
      goto 1000
   endif
! a TPfun can have different number of ranges, must be allocated
   displace=4+nwch(16)
   allocate(tpfuns(jfun)%limits(nr))
   allocate(tpfuns(jfun)%funlinks(nr))
   call loadrn(nr,iws(lfun+displace),tpfuns(jfun)%limits)
   call loadr(lfun+displace+nr*nwpr,iws,tpfuns(jfun)%hightlimit)
!   write(*,*)'3Z high T',tpfuns(jfun)%hightlimit
! the expressions are linked from here, one per range
   displace=displace+(1+nr)*nwpr
   loklexpr=lfun+displace-1
! extract  the expressions
   do kx=1,nr
      lexpr=iws(loklexpr+kx)
      nc=iws(lexpr)
!      write(*,*)'3Z coeffs 1',kx,lfun,lexpr,nr,nc
      exprot=>tpfuns(jfun)%funlinks(kx)
      exprot%noofcoeffs=nc
!      if(nc.gt.20) stop
      allocate(exprot%tpow(nc))
      allocate(exprot%ppow(nc))
      allocate(exprot%wpow(nc))
      allocate(exprot%plevel(nc))
      allocate(exprot%link(nc))
      allocate(exprot%coeffs(nc))
      mmz=lexpr
!      write(*,*)'3Z coeffs 2',nc,iws(nc),mmz
      do i=1,nc
         exprot%link(i)=iws(mmz+1)
         exprot%tpow(i)=iws(mmz+2)
         exprot%ppow(i)=iws(mmz+3)
         exprot%wpow(i)=iws(mmz+4)
         exprot%plevel(i)=iws(mmz+5)
         call loadr(mmz+6,iws,exprot%coeffs(i))
         mmz=mmz+5+nwpr
      enddo
   enddo
!   if(jfun.gt.0) then
!      read(lut)tpfuns(jfun)%hightlimit
!   else
!      read(lut)dummy
!   endif
1000 continue
   return
 end subroutine read0tpfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine makeoptvname
!\begin{verbatim}
 subroutine makeoptvname(name,indx)
    implicit none
    character name*(*)
    integer indx
!\end{verbatim} %+
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
!1000 continue
    return
  end subroutine makeoptvname

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine findtpused
!\begin{verbatim}
 subroutine findtpused(lfun,string)
! this routine finds which other TPFUNS (including parameters) that
! use the TPFUN lfun.  It is used when listing optimizing coefficients
   implicit none
   integer lfun
   character string*(*)
!\end{verbatim} %+
   integer jp,kfun,nr,nc,j1
   type(tpfun_expression), pointer :: exprot
   string=' '
   jp=1
   loop1: do kfun=1,freetpfun-1
      if(kfun.eq.lfun) cycle
      loop2: do nr=1,tpfuns(kfun)%noofranges
         exprot=>tpfuns(kfun)%funlinks(nr)
         if(.not.associated(exprot)) cycle loop2
         nc=exprot%noofcoeffs
         loop3: do j1=1,nc
            if(exprot%link(j1).eq.lfun) then
!               write(*,*)'3Z found: ',trim(tpfuns(kfun)%symbol),kfun
               string(jp:)=tpfuns(kfun)%symbol
               jp=len_trim(string)+2
               cycle loop1
            endif
         enddo loop3
      enddo loop2
   enddo loop1
!   write(*,*)'3Z where: ',trim(string)
!1000 continue
   return
 end subroutine findtpused

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine list_tpfun_details
!\begin{verbatim}
 subroutine list_tpfun_details(lfun)
! listing the internal datastructure of all tpfuns
! converts all TP functions to arrays of coefficients with powers of T
   implicit none
   integer lfun
!\end{verbatim}
   integer j1,j2,j3,nc
   TYPE(tpfun_expression), pointer :: exprot
   if(lfun.lt.0) then
! list all ...
      continue
   elseif(lfun.ge.freetpfun) then
      write(*,*)'No such function'
   else
      exprot=>tpfuns(lfun)%funlinks(1)
      nc=exprot%noofcoeffs
      write(*,100)tpfuns(lfun)%symbol,tpfuns(lfun)%noofranges,nc,&
           firsteq%eq_tpres(lfun)%results(1)
100   format('Name: ',a,2i5,(1pe12.4)/&
           '    term  coefficent    tpow  ppow  wpow plevel  link')
      do j1=1,nc
         write(*,110)j1,exprot%coeffs(j1),exprot%tpow(j1),exprot%ppow(j1),&
              exprot%wpow(j1),exprot%plevel(j1),exprot%link(j1)
110      format('Term: ',i2,1pe12.4,2x,5i6)
      enddo
   endif
!1000 continue
   return
 end subroutine list_tpfun_details

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!
! Below are a couple of routines to generate SOLGASMIX DAT files
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
 
!\addtotable subroutine tpfun2coef
!\begin{verbatim}
 subroutine tpfun2coef(ctpf,ntpf,npows,text)
! called by saveadatformat in gtp3C to generate SOLGASMIX DAT files
! converts all TP functions to arrays of coefficients with powers of T
   implicit none
   integer ntpf,npows
   type(gtp_tpfun2dat) :: ctpf(*)
   character text*(*)
!\end{verbatim} %+
! powers are 0  1  100    2  3  -1 ; 7  -9  -2  -3  extra
!                  Tln(T)            these on extra line
   integer, parameter :: maxnc=15
   integer i1,i2,i3,usedpow(maxnc)
   character buffer*80
   logical done
!   write(*,*)'In tpfun2coef with ',ntpf,' ctpf records allocated'
   do i1=1,ntpf
      ctpf(i1)%nranges=-1
   enddo
! this loop may have to be done several times as functions calling functions
   done=.false.
   do while(.not.done)
      done=.true.
! skip the first two functions ... R and RTLNP
      do i1=3,ntpf
! done is set false if the function ctpf(i1) is not converted
         call tpf2c(ctpf,i1,done)
         if(gx%bmperr.ne.0) goto 1000
      enddo
   enddo
! here all TP functions are converted to coefficients
!   i1=19
!   call tpwrite('ee',i1,ctpf(i1)%nranges,ctpf(i1)%cfun)
! extract the powers used
   npows=0
   do i1=3,ntpf
      do i2=1,ctpf(i1)%nranges
         call sortcoeffs(maxnc,i1,ctpf(i1)%cfun%coefs(1,i2),&
              ctpf(i1)%cfun%tpows(1,i2))
         call checkpowers(maxnc,i1,ctpf(i1)%cfun%tpows(1,i2),npows,usedpow)
      enddo
   enddo
   buffer=' '
!   write(*,12)npows,(usedpow(i1),i1=1,npows)
!   write(buffer,11)npows,(usedpow(i1),i1=1,npows)
!11 format(12i5)
!12 format('3Z power2: ',i3,12i4)
   text=buffer
!   do i1=3,ntpf
!      call tpwrite('ee',i1,ctpf(i1)%nranges,ctpf(i1)%cfun)
!      write(*,*)'Sorting function/range: ',i1,ctpf(i1)%nranges
!      write(*,699)i1,ctpf(i1)%nranges
!      do i2=1,ctpf(i1)%nranges
!         write(*,700)ctpf(i1)%cfun%tbreaks(i2),&
!              (ctpf(i1)%cfun%coefs(i3,i2),i3=1,npows)
!      enddo
!   enddo
!800 format(a,2i3,3(1pe12.4,i5))
!699 format('Function/parameter and ranges: ',2i4)
!700 format(F11.4,4x,4(1x,G14.8)/5(1x,G14.8))
1000 continue
   return
 end subroutine tpfun2coef

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine list_tpascoef
!\begin{verbatim} %-
 subroutine list_tpascoef(lut,text,paratyp,i1,npows,factor,ctpf)
! writes a parameter in DAT format
! text contains the stoichiometries written with the format 1x,F11.6
! it can be very long if there are many coefficients.
   implicit none
   integer lut,i1,npows,paratyp
   character text*(*)
! this is a factor that may be multiplied with all coefficients for
! phases like sigma which has only a disordered part. Also ionic liquid
   double precision factor
   type(gtp_tpfun2dat) :: ctpf(*)
!\end{verbatim} %+
   integer i2,i3,ip,kk,mm
   ip=len_trim(text)
! this is the endmember stoichiometry, 12 characters per value, 6x12=72
!   write(*,*)'3Z len_trim(text): ',ip
   if(ip.gt.72) then
      write(lut,698)paratyp,ctpf(i1)%nranges,text(1:72)
      i2=73
      do while(i2.lt.ip)
         write(lut,699)trim(text(i2:i2+71))
         i2=i2+72
      enddo
   else
      write(lut,698)paratyp,ctpf(i1)%nranges,trim(text)
   endif
!698 format(i4,i3,a)
! According to Ted
698 format(i2,i3,1x,a)
699 format(a)
   do i2=1,ctpf(i1)%nranges
      if(ctpf(i1)%cfun%coefs(7,i2).eq.zero .and. &
           ctpf(i1)%cfun%coefs(8,i2).eq.zero .and. &
           ctpf(i1)%cfun%coefs(9,i2).eq.zero .and. &
           ctpf(i1)%cfun%coefs(10,i2).eq.zero .and. &
           ctpf(i1)%cfun%coefs(11,i2).eq.zero) then
         write(lut,700)ctpf(i1)%cfun%tbreaks(i2),&
              (factor*ctpf(i1)%cfun%coefs(i3,i2),i3=1,6)
      else
! There are some special powers, write only non-zero coefficients
!         write(lut,705)ctpf(i1)%cfun%tbreaks(i2),&
!              (ctpf(i1)%cfun%coefs(i3,i2),i3=1,6),&
!              (ctpf(i1)%cfun%coefs(i3,i2),&
!              ctpf(i1)%cfun%tpows(i3,i2),i3=7,npows)
         write(lut,710)ctpf(i1)%cfun%tbreaks(i2),&
              (factor*ctpf(i1)%cfun%coefs(i3,i2),i3=1,6)
         mm=0
! The 6 first powers are the default 0 1 100 2 3 -1
! Possible extra powers are 7 -9 -2 unknown1 unknown2         
! uknown can be -3, 4, 5, -8 (for sqrt(T),
         do kk=7,npows
            if(ctpf(i1)%cfun%coefs(kk,i2).ne.zero) mm=mm+1
         enddo
!         write(*,719)'3Z powers: ',npows,mm,&
!           (ctpf(i1)%cfun%coefs(i3,i2),ctpf(i1)%cfun%tpows(i3,i2),i3=7,npows)
! UNFINISHED
         write(lut,730,advance='no')mm
         do kk=7,npows
            if(ctpf(i1)%cfun%coefs(kk,i2).ne.zero) then
               if(mm.eq.1) then
                  if(ctpf(i1)%cfun%tpows(kk,i2).eq.-8) then
! this is the square root of T
                     write(lut,733)factor*ctpf(i1)%cfun%coefs(kk,i2)
                  else
                     write(lut,731)factor*ctpf(i1)%cfun%coefs(kk,i2),&
                          ctpf(i1)%cfun%tpows(kk,i2)
                  endif
                  mm=mm-1
               elseif(mm.lt.0) then
                  write(*,*)'3Z wrong number of coefficeints!!!'
               else
                  if(ctpf(i1)%cfun%tpows(kk,i2).eq.-8) then
! this is the square root of T
                     write(lut,733)factor*ctpf(i1)%cfun%coefs(kk,i2)
                  else
                     write(lut,731,advance='no')&
                          factor*ctpf(i1)%cfun%coefs(kk,i2),&
                          ctpf(i1)%cfun%tpows(kk,i2)
                  endif
                  mm=mm-1
               endif
            endif
         enddo
      endif
   enddo
! according to Ted
700 format(1x,F11.4,6(1x,G15.8)/' 1 0.00000000       0.00')
!705 format(1x,F11.4,6(1x,G15.8)/' 3 ',3(1x,G15.8,i3,'.00'))
710 format(1x,F11.4,6(1x,G15.8))
!719 format(a,2i3,4(1x,G10.2,1x,i3,'.00'))
!721 format(1x,i3,4(1x,G15.8,1x,F5.2))
730 format(i3)
731 format(1x,G15.8,1x,i3,'.00')   
733 format(1x,G15.8,1x,'  0.50')
!1000 continue
   return
 end subroutine list_tpascoef

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine tpf2c
!\begin{verbatim} %-
 subroutine tpf2c(ctpf,lfun,done)
! convert TPfun lfun to an array of coefficients with powers of T
! if this TP function already converted just return
! if this TP function calls another TP function not converted return error
!
   implicit none
   integer lfun
   logical done
   type(gtp_tpfun2dat) :: ctpf(*)
!\end{verbatim} %+
   integer i1,i2,i3,nrange,nc,funref
   type(tpfun_root), pointer :: tpfroot
   type(tpfun_expression), pointer :: tpfexpr
! return if already converted
!   write(*,*)'3Z in tpf2c ',lfun,tpfuns(lfun)%noofranges
   if(ctpf(lfun)%nranges.ge.0) goto 1000
! This function not converted, check if it reference an unconverted TPfunction
   tpfroot=>tpfuns(lfun)
   if(btest(tpfuns(lfun)%status,TPCONST)) then
      write(*,*)'3Z this function is a constant, have to think about',lfun
      stop 18
   endif
   nrange=tpfroot%noofranges
   do i1=1,nrange
      tpfexpr=>tpfroot%funlinks(i1)
      nc=tpfexpr%noofcoeffs
! skip the first two predefined functions, R and RTLNP
      do i2=1,nc
         funref=tpfexpr%link(i2)
         if(funref.eq.1) then
! this is a constant R, multiply the coefficient with 8.31451 and set link=0
!            write(*,*)'3Z Replacing R with its value in function ',lfun
            tpfexpr%coeffs(i2)=8.31451*tpfexpr%coeffs(i2)
            tpfexpr%link(i2)=0
         elseif(funref.eq.2) then
!            write(*,*)'3Z Deleting use of RTLNP for gas in function ',lfun
            tpfexpr%link(i2)=0
            tpfexpr%coeffs(i2)=zero
         elseif(funref.gt.0) then
            if(ctpf(funref)%nranges.lt.0) then
! this function has a reference to an unconverted TPfunction
!               write(*,*)'3Z TPfun ',lfun,' reference ',funref,&
!                    ctpf(funref)%nranges
               done=.false.
!               write(*,*)'Skipping for the moment ',lfun,funref
               goto 1000
            endif
         endif
      enddo
   enddo
! convert the TPfun "lfun" to coefficents and powers
   call tpf2cx(ctpf,lfun,nrange,ctpf(lfun)%cfun)
!   write(*,*)'3Z tp2c: ',tpfuns(lfun)%symbol,nrange
!   call tpwrite('z2',lfun,nrange,ctpf(lfun)%cfun)
!   do i1=1,nrange
!      write(*,200)(ctpf(lfun)%cfun%coefs(i2,i1),i2=1,6)
!   enddo
!200 format(10F12.3)
1000 continue
   return
 end subroutine tpf2c

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

! each term has a coefficent and an array of integers
! tpow is power of T
! ppow is power of P
! wpow is power of linked symbol, the link is in link
! plevel is level of parenthesis ??
! link is link to another function if >0 or a unary function if <0
!      accept only -2 which is taken as LN(T)
!\addtotable subroutine tpf2cx
!\begin{verbatim} %-
 subroutine tpf2cx(ctpf,lfun,nrange,cfun1)
! convert TPfun lfun to an array of coefficients with powers of T
! if this TP function already converted just return
! if this TP function calls another TP function not converted return error
   implicit none
   integer lfun,nrange
   type(gtp_tpfun_as_coeff) :: cfun1
   type(gtp_tpfun2dat) :: ctpf(*)
!   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
! max no of coefficent, max no of ranges ...
   integer, parameter :: maxnc=15,maxnr=20
   integer i1a,i1b,i2,i3,i4,nc1,funref,iadd,nrangeb,ncc,caddnr(maxnr),nrr,jadd
   integer caddid(maxnr),krange,klink
   type(tpfun_root), pointer :: tpfroot
   type(tpfun_expression), pointer :: tpfexpr
   type(gtp_tpfun_as_coeff), dimension(:), allocatable :: cadd
   double precision ccc
   logical skipnext,isqrt
   integer noofcadd,fsqrt
!----------
!  TYPE gtp_tpfun_as_coeff
! record for a TPFUN converted to coefficents without any references to other
! functions.  Note ranges may change when adding functions!!
!     double precision, dimension(:), allocatable :: tbreaks
!     double precision, dimension(:,:), allocatable :: coefs
!     integer, dimension(:,:), allocatable :: tpows
!  end type gtp_tpfun_as_coeff
!----------
! Now convert !!
! allocate a record for all ranges and coefficients
   iadd=0
   tpfroot=>tpfuns(lfun)
   if(nrange.gt.maxnr) then
      write(*,*)'3Z too many T ranges!',nrange
      stop 13
   endif
!   write(*,*)'3Z converting ',lfun
! allocate and zero cfun1 data
   allocate(cfun1%tbreaks(maxnr))
   allocate(cfun1%coefs(maxnc,maxnr))
   allocate(cfun1%tpows(maxnc,maxnr))
   cfun1%tbreaks=zero
   cfun1%coefs=zero
   cfun1%tpows=-100
! T ranges, low T limit ignored.  NOTE the number of ranges may change!
   do i1a=2,nrange
      cfun1%tbreaks(i1a-1)=tpfroot%limits(i1a)
   enddo
   cfun1%tbreaks(nrange)=tpfroot%hightlimit
   nrangeb=nrange
!   write(*,2)tpfuns(lfun)%symbol
!2  format(/'Entering tpf2cx ---------------------------------: ',a)
! functions
! NOTE nrange may change below if referenced functions have a smaller range
! not so good to have a loop ...
! i1a is range index for original TPfun
! i1b is range index for TPfun converted to coefficients as referenced
!                        functions may have shorter range
!   do i1b=1,nrange
!      tpfexpr=>tpfroot%funlinks(i1b)
!      nc1=tpfexpr%noofcoeffs
!      do i1a=1,nc1
!         write(*,30)i1a,tpfexpr%coeffs(i1a),tpfexpr%tpow(i1a),&
!              tpfexpr%ppow(i1a),tpfexpr%wpow(i1a),tpfexpr%plevel(i1a),&
!              tpfexpr%link(i1a)
!      enddo
!   enddo
!30 format('3Z term: ',i3,1pe12.4,5i7)
   i1a=0
   i1b=0
   jadd=0
   noofcadd=0
   krange=1
   isqrt=.false.
100 continue
   i1a=i1a+1
   i1b=i1b+1
   if(i1a.gt.nrange) goto 700
      jadd=jadd+i1a
      tpfexpr=>tpfroot%funlinks(i1a)
      nc1=tpfexpr%noofcoeffs
!      write(*,107)'TPfun ranges: ',tpfuns(lfun)%symbol,i1a,i1b,nrangeb,jadd
!107   format(a,a,5i4)
      skipnext=.false.
! maybe needed? YES!
      iadd=0
      trange: do i2=1,nc1
         cfun1%coefs(i2,i1b)=tpfexpr%coeffs(i2)
         cfun1%tpows(i2,i1b)=tpfexpr%tpow(i2)
! assume link to unary LN function just means LN(T)
         funref=tpfexpr%link(i2)
!         write(*,113)'3Z cc: ',i1a,i1b,i2,funref,&
!              cfun1%coefs(i2,i1b),cfun1%tpows(i2,i1b)
!113      format(a,4i4,E20.8,i5)
         if(skipnext) then
! skip this term as it should just contain the ln(T)
            if(tpfexpr%plevel(i2).ne.1 .or. tpfexpr%link(i2).ne.0 &
                 .or. tpfexpr%tpow(i2).ne.1) then
!               write(*,*)'3Z WARNING check if TPFUN error in: ',&
!                    trim(tpfroot%symbol),lfun
!               gx%bmperr=4393; goto 1000
            endif
            cfun1%coefs(i2,i1b)=zero
            cfun1%tpows(i2,i1b)=-100
            skipnext=.false.
            if(isqrt) then
! fixing a bug for T**2 as EXP(0.5*LN(T))
               isqrt=.false.
               cycle trange
            endif
         endif
         funrefif: if(funref.lt.0) then
! this is assumed to be a link to LN(T) or SQRT
            if(funref.eq.-3) then
! this is to handle sqrt(t) which in a TDB file is EXP(0.5LN(T))
! check that link to function is SQRT
               fsqrt=tpfexpr%link(i2+1)
!               write(*,114)'Found SQRTT?: ',lfun,tpfuns(lfun)%symbol,&
!                    fsqrt,tpfuns(fsqrt)%symbol
!114            format(a,i5,2x,a,i5,2x,a)
! set T power to -8; this will be converted to 0.5 when writing !!!
               cfun1%tpows(i2,i1b)=-8
               write(*,*)'3Z sqrt coeff: ',cfun1%coefs(i2,i1b)
               isqrt=.true.
               skipnext=.true.
            elseif(funref.ne.-2) then
! this is an unknown type of funref link (-2 means LN(T))
               write(*,*)'3Z TPFUN with other unary function than LN(T): ',&
                    trim(tpfroot%symbol),lfun,funref
               gx%bmperr=4393; goto 1000
            elseif(tpfexpr%tpow(i2).eq.1) then
! NOTE not elseif(funref ... just extract the power of T, could be 0?
! Tln(T) will have tpows = 100, we skip the next term with the T
               cfun1%tpows(i2,i1b)=tpfexpr%tpow(i2)+99
               skipnext=.true.
!            else
! This could be a LN(T) term?
!               write(*,'(a,a,5i5)')'3Z TPFUN with just LN(T)? ',&
!                    trim(tpfroot%symbol),lfun,i2,funref,tpfexpr%tpow(i2)
!               gx%bmperr=4393; goto 1000
            endif
         elseif(funref.gt.0) then
            funrefranges: if(ctpf(funref)%nranges.gt.0) then
! this range has a reference to a converted TPfunction,
! store this separately, possibly multiplied with coefficent and T powers
! and link all such functions to be added using cfun1%nextcrec
! examples:  +22*GHSERCR, ff*exp(qq*irt) ... the latter will not work ...
               ccc=tpfexpr%coeffs(i2)
!               write(*,32)'3Z link from: ',tpfuns(lfun)%symbol,&
!                    i2,funref,ccc,trim(tpfuns(funref)%symbol)
!32             format(a,a,2i4,F6.2,' to ',a)
!               write(*,333)'3Z term, link, factor: ',i2,funref,ccc,&
! IMPORTRANT funref is also index in tpfuns!!!
!                    trim(tpfroot%symbol),trim(tpfuns(funref)%symbol),&
!                    ctpf(funref)%nranges,&
!                    tpfuns(funref)%noofranges,size(tpfuns(funref)%limits)
!333            format(a,2i4,1pe11.2,2x,a,2x,a,5i5)
! only allow a constant coefficent, no T or P powers, no unary function ...
               if(tpfexpr%tpow(i2).ne.0 .or. tpfexpr%ppow(i2).ne.0 .or. &
                    tpfexpr%wpow(i2).ne.0 .and. &
                    (tpfexpr%plevel(i2).ne.0 .or. tpfexpr%plevel(i2).ne.1)) then
! Above the function SQRT which has tpfexpr%plevel(i2)=1 is accepted ...
                  if(tpfuns(funref)%noofranges.eq.1) then
! Now check if funref is just a constant, then multiply ccc with that!
!                     write(*,334)'3Z trying to handle MEV factor ... ',&
!                          trim(tpfuns(funref)%symbol),&
!                          tpfuns(funref)%funlinks(1)%noofcoeffs,&
!                          tpfuns(funref)%funlinks(1)%coeffs(1)
!334                  format(a,a,i2,1pe11.2)
! WOW wpow-1000 is link to another function!
                     klink=tpfexpr%wpow(i2)-1000
                     if(tpfuns(funref)%funlinks(1)%noofcoeffs.eq.1 .and.&
                          tpfuns(klink)%funlinks(1)%noofcoeffs.eq.1) then
                        ccc=ccc*tpfuns(funref)%funlinks(1)%coeffs(1)*&
                             tpfuns(klink)%funlinks(1)%coeffs(1)
!                        write(*,335)'3Z Wow! ',trim(tpfuns(funref)%symbol),&
!                             trim(tpfuns(lfun)%symbol),i2,nc1,ccc,&
!                             tpfuns(funref)%funlinks(1)%coeffs(1),&
!                             klink,trim(tpfuns(klink)%symbol),&
!                             tpfuns(klink)%funlinks(1)%coeffs(1)
!335                     format(a,2x,a,2x,a,2i3,2(1pe12.4),i3,2x,a,1pe12.4)
!                        cfun1%coefs(i2,i1b)=ccc
                        cfun1%coefs(i2,i1b)=ccc
                        exit funrefif
                     endif
                  endif
! else give up
                  write(*,116)'3Z Too complicated function: ',&
                       trim(tpfroot%symbol),tpfexpr%tpow(i2),&
                       tpfexpr%ppow(i2),tpfexpr%wpow(i2),tpfexpr%plevel(i2),&
                       funref,trim(tpfuns(funref)%symbol)
116               format(a,a,5i5,2x,a)
                  gx%bmperr=4399; goto 1000
               endif
! this term should be ignored as it replaced by the function
               cfun1%coefs(i2,i1b)=zero
               cfun1%tpows(i2,i1b)=-100
! we must create a new coefficient array with the funref coefficents
! multiplied with the current coef within the current T-range
! It may be necessary to increase the number of T-ranges
               if(.not.allocated(cadd)) then
! we have more than 6 functions added in soma cases ...
                  allocate(cadd(10))
                  noofcadd=noofcadd+1
!                  write(*,*)'Allocating cadd ',i2,noofcadd
! ??                  iadd=0
                  caddnr=0
               endif
! call a new function to add the coefficents of funref
! to ctpf(
               iadd=iadd+1
               if(iadd.gt.7) then
                  write(*,*)'3Z many added functions in: ',&
                       trim(tpfuns(lfun)%symbol),': ',iadd,funref
               endif
               cadd(iadd)=ctpf(funref)%cfun
               caddnr(iadd)=ctpf(funref)%nranges
               caddid(iadd)=funref
!               write(*,800)'3Z aa: ',funref,iadd,(cadd(iadd)%coefs(i3,1),&
!                    cadd(iadd)%tpows(i3,1),i3=1,3)
! multiply all terms in funref with the coefficient of this term
! within the current T-range               
! It may be necessary to increase the T-ranges of ctpf
!               write(*,*)'3Z addranges: ',ctpf(funref)%nranges,ccc
               do i3=1,maxnc
                  do i4=1,ctpf(funref)%nranges
                     cadd(iadd)%coefs(i3,i4)=ccc*cadd(iadd)%coefs(i3,i4)
                  enddo
               enddo
!               write(*,800)'3Z bb: ',funref,iadd,(cadd(iadd)%coefs(i3,1),&
!                    cadd(iadd)%tpows(i3,1),i3=1,3)
            else
! what about funref with no ranges?
               write(*,*)'3Z funref has no ranges? ',&
                    trim(tpfuns(funref)%symbol),ctpf(funref)%nranges
               gx%bmperr=4399; goto 1000
            endif funrefranges
!         else
! when funref=0 it is OK to do nothing !!
         endif funrefif
!800      format(a,2i3,3(1pe12.4,i5))
! we have gone through all terms for the TPfun for this range
      enddo trange
! Check if there were function links in this range
      if(iadd.gt.0) then
! If iadd>1 we must first add together all the different functions referenced
! and possibly split the T range if these function have a different ranges
         ncc=3
!         write(*,*)'3Z adjusting ranges?',iadd
! This loop only if there are two or more function references within a range
         do i3=iadd,2,-1
            nrr=caddnr(i3)
!            write(*,16)'3Z there are coefficients to add!!',i3,iadd,nrr
!16          format(a,6i4)
! add terms and adjust all ranges in cadd(i3-1)
            call adjustranges(nrr,cadd(i3-1),caddnr(i3),cadd(i3))
            if(gx%bmperr.ne.0) then
               write(*,*)'Error occured adding: ',caddid(i3-1),caddid(i3)
               goto 1000
            endif
! note nrr may be updatad
            caddnr(i3-1)=nrr
! we may have more functions to add ...
         enddo
!         call tpwrite('++',0,caddnr(1),cadd(1))
! we have now added all links, now add the sum of all cadd to cfun1 range i1a
! adjust1ranges creates breakpoints only in the current range, i1a, of cfun1
! and adds the coefficients from cadd(1) to this
!         write(*,*)'3Z calling adjust1: ',i1a,jadd,i1b,nrangeb
!         call tpwrite('>1',lfun,nrangeb,ctpf(lfun)%cfun)
!         call tpwrite('>2',0,caddnr(1),cadd(1))
         call adjust1range(lfun,jadd,nrangeb,krange,cfun1,caddnr(1),cadd(1))
         krange=krange+1
! if additional ranges needed nrangeb changed 
! increment i1b but not i1a and nrange
! why -1 ??
         jadd=nrangeb-i1a-1
         i1b=i1b+nrangeb-2
!         write(*,*)'3Z after adjust1x: ',i1a,jadd,i1b,nrangeb
!         call tpwrite('<<',lfun,nrangeb,ctpf(lfun)%cfun)
!         write(*,*)'deallocating cadd'
         deallocate(cadd)
      endif
      goto 100
! we have gone through all ranges
700 continue
!   write(*,*)'3Z 700: ',i1a,nrange,nrangeb
   ctpf(lfun)%nranges=nrangeb
!   write(*,*)'3Z converted function with ranges: ',lfun,nrangeb
! Listing the final function
!   call tpwrite('z1',lfun,nrange,ctpf(lfun)%cfun)
1000 continue
   return
 end subroutine tpf2cx

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

 subroutine tpwrite(c2,lfun,nrange,cfun)
! temporary debug output
! ************************* ATTENTION
! IF I DO NOT CALL THIS THERE ARE BUGS !!
   IMPLICIT none
   integer nrange,lfun
   character c2*2
   type(gtp_tpfun_as_coeff) :: cfun
! 
   integer i1,i2
   do i1=1,nrange
!      write(*,800)c2,lfun,i1,nrange,cfun%tbreaks(i1),&
!           (cfun%coefs(i2,i1),cfun%tpows(i2,i1),i2=1,10)
   enddo
!800 format('3Z ',a,': ',i3,2i2,F9.2,3(1pe13.5,i5)/&
!         (23x,e13.5,i5,e13.5,i5,e13.5,i5))
!   write(*,*)'3Z end of function'
   return
 end subroutine tpwrite

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

 subroutine tpmult(lfun,mfun,ccc,ctpf)
! multiples all terms in cfpf(lfun) with the factor ccc and returns that
! in ctpf(mfun)
! cfun is not changed
   implicit none
   integer lfun,mfun
   double precision ccc
   type(gtp_tpfun2dat) :: ctpf(*)
!   type(gtp_tpfun_as_coeff) :: cfun
! 
   integer, parameter :: maxnc=15,maxnr=20
   integer i1,i2
   ctpf(mfun)%nranges=ctpf(lfun)%nranges
   if(.not.allocated(ctpf(mfun)%cfun%tbreaks)) then
!      write(*,*)'3Z allocating mfun'
      allocate(ctpf(mfun)%cfun%tbreaks(maxnr))
      allocate(ctpf(mfun)%cfun%coefs(maxnc,maxnr))
      allocate(ctpf(mfun)%cfun%tpows(maxnc,maxnr))
   endif
   ctpf(mfun)%cfun%tbreaks=zero
   ctpf(mfun)%cfun%coefs=zero
   ctpf(mfun)%cfun%tpows=0
   do i1=1,ctpf(mfun)%nranges
      ctpf(mfun)%cfun%tbreaks(i1)=ctpf(lfun)%cfun%tbreaks(i1)
      do i2=1,maxnc
         ctpf(mfun)%cfun%coefs(i2,i1)=ccc*ctpf(lfun)%cfun%coefs(i2,i1)
         ctpf(mfun)%cfun%tpows(i2,i1)=ctpf(lfun)%cfun%tpows(i2,i1)
      enddo
   enddo
   return
 end subroutine tpmult

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine adjust1range
!\begin{verbatim} %-
 subroutine adjust1range(lfun,nr1,nrange,krange,ctp1,nr2,ctp2)
! check if ctp1 range nr1 must be split in more ranges due to tbreaks in ctp2
! nrange is the total number of ranges of ctp1
! There are 10 ranges allocated for all, nr1 and nr2 are the used ranges
   implicit none
   integer lfun,nr1,nr2,krange,nrange
   type(gtp_tpfun_as_coeff) :: ctp1,ctp2
!\end{verbatim} %+
!  TYPE gtp_tpfun_as_coeff
! record for a TPFUN converted to coefficents without references to other funs
! Note ranges may increase when adding functions!!
!     double precision, dimension(:), allocatable :: tbreaks
!     double precision, dimension(:,:), allocatable :: coefs
!     integer, dimension(:,:), allocatable :: tpows
!  end type gtp_tpfun_as_coeff
   double precision, parameter :: tenth=1.0D-1
   integer, parameter :: maxnc=15
   integer i1,i2,i3,k1,nr3,nr0,j2,j3,mrange
   double precision tlow1,thigh1,tlow2,thigh2
   logical nosplit
   type(gtp_tpfun_as_coeff) :: ctp3
!
   nr0=nr1
   if(nr1.eq.1) then
      tlow1=298.15
      thigh1=ctp1%tbreaks(nr1)
      j2=1
   else
      tlow1=ctp1%tbreaks(nr1-1)
      thigh1=ctp1%tbreaks(nr1)
      j2=nr1-1
   endif
!   write(*,7)'adjust1range: ',nrange,nr1,nr2,j2,&
!        ctp1%tbreaks(nr1),ctp2%tbreaks(nr2)
!7  format(a,4i4,2F10.2)
   allocate(ctp3%tbreaks(1))
   allocate(ctp3%coefs(maxnc,1))
   allocate(ctp3%tpows(maxnc,1))
   ctp3%tbreaks=zero
   nosplit=.true.
! search ctp2 for tbreaks in the range tlow1 to thigh1
   i2=1
!   mrange=nrange-1
!   write(*,16)'In adjust1range for range: ',nr1,i2,krange,&
!        tlow1,thigh1,ctp2%tbreaks(i2)
!16 format(/a,3i3,3F10.2)
!100 continue
   split: do while(i2.lt.nr2)
! fine-tuning needed when breakpoints identical in parameter and GHSERxx
!      write(*,*)'in do while: ',i2,nr2,ctp2%tbreaks(i2),tlow1
      if(ctp2%tbreaks(i2)-tlow1.gt.tenth) then
!         write(*,16)'3Z check breakpoint ',nrange,i2,krange,&
!              ctp2%tbreaks(i2),thigh1
! fine-tuning needed when breakpoints identical in parameter and GHSERxx
         if(abs(ctp2%tbreaks(i2)-thigh1).lt.tenth) then
! breakpoints are identical
            mrange=nrange-1
!            write(*,*)'Identical breakpoints',ctp2%tbreaks(i2),thigh1
            goto 800
         elseif(ctp2%tbreaks(i2)-thigh1.lt.-tenth) then
! fine-tuning needed when breakpoints identical in parameter and GHSERxx
! there is a breakpoint in ctp2 between tlow1 and thigh1
! we must add one range above nr1, shift the coefficients in higher ranges up 
!            write(*,16)'3Z new breakpoint ',nrange,j2,0,ctp2%tbreaks(i2),thigh1
!            call tpwrite('--',0,nrange,ctp1)
            do k1=nrange,j2,-1
               ctp1%tbreaks(k1+1)=ctp1%tbreaks(k1)
               do i3=1,maxnc
! copy the coefficients to the new range
                  ctp1%coefs(i3,k1+1)=ctp1%coefs(i3,k1)
                  ctp1%tpows(i3,k1+1)=ctp1%tpows(i3,k1)
               enddo
            enddo
            nrange=nrange+1
!            call tpwrite('up',0,nrange,ctp1)
! now add coeffs from ctp1 range k1 and ctp2 in range i2 to ctp3 range 1
! then replace range k1 in ctp1 by range 1 of ctp3
            ctp3%tpows=-100
            ctp3%coefs=zero
! the range in ctp1 that should be added to is j2+1 ??
!            j2=j2+1
!            write(*,*)'3Z add7: ',nr1,i2,nrange,j2
!            call tpwrite('v1',0,nrange,ctp1)
!            call tpwrite('v2',0,nr2,ctp2)
            call add1tpcoeffs(j2,ctp1,i2,ctp2,1,ctp3)
!            call tpwrite('v3',0,1,ctp3)
            do j3=1,maxnc
               ctp1%coefs(j3,j2)=ctp3%coefs(j3,1)
               ctp1%tpows(j3,j2)=ctp3%tpows(j3,1)
            enddo
! NEW: added a range to ctp1 !!!
            krange=krange+1
            tlow1=min(ctp2%tbreaks(i2),thigh1)
            ctp1%tbreaks(j2)=tlow1
!            call tpwrite('q1',0,nrange,ctp1)
! we have added one range to ctp1
            nosplit=.false.
!            nrange=nrange+1
!            write(*,*)'adjust1range 6:',nrange,j2,ctp1%tbreaks(j2)
         else
            mrange=nrange-1
            goto 800
         endif
      else
         continue
!         write(*,731)'adjust1range 7:',i2,nr2,ctp2%tbreaks(i2),tlow1
!731      format(a,2i3,2F10.2)
      endif
      i2=i2+1
      j2=j2+1
   enddo split
! flyttat till efter label 800
!799 continue
   mrange=nrange
800 continue
!   write(*,*)'Why??',nrange,mrange,nr1,nr2,nosplit
   call tpwrite('w0',0,nrange,ctp1)
! just add the terms (for the last range)
!   ctp3%tpows=-100
!   ctp3%coefs=zero
   if(nosplit) then
! we have not split the range, store cp3 in range nr1
      mrange=nr1
   endif
!   write(*,900)'3Z add8: ',nrange,i2,mrange,krange,&
!        ctp1%tbreaks(nrange),ctp2%tbreaks(i2)
!900 format(/a,4i3,2F10.2)
! these output w1..c4 are important for debugging
   call tpwrite('w1',0,nrange,ctp1)
   call tpwrite('w2',0,nr2,ctp2)
!   call add1tpcoeffs(mrange,ctp1,i2,ctp2,1,ctp3)
   call add1tpcoeffs(nrange,ctp1,i2,ctp2,1,ctp3)
   call tpwrite('w3',0,1,ctp3)
! skipping this loop
   if(krange.ne.nrange) then
      write(*,*)' *** Check function: ',tpfuns(lfun)%symbol,nrange,krange
   endif
! We need to know which range in ctp1 we should store ctp3 .... krange!!
! to handle problems with G(LIQ,U+4:O-2) parameter
   do j3=1,maxnc
      ctp1%coefs(j3,krange)=ctp3%coefs(j3,1)
      ctp1%tpows(j3,krange)=ctp3%tpows(j3,1)
   enddo
   call tpwrite('w4',0,nrange,ctp1)
!1000 continue
!   if(nr3.gt.nr0) then
!      write(*,*)'3Z inserted ',nrange-nr0,' ranges'
!   endif
   return
 end subroutine adjust1range

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine adjustranges
!\begin{verbatim} %-
 subroutine adjustranges(nr1,ctp1,nr2,ctp2)
! add ctp2 to ctp1 which to have the same ranges and breakpoints
! nr1 and nr2 give the number of T-ranges and breakpoints
! add coefficients of ctp1 and ctp2 for each range
! NOTE these already multiplied with the coefficents!!
! There are 10 coefficients allocated for all functions.
! the added function is returned as ctp1
   implicit none
   integer nr1,nr2
   type(gtp_tpfun_as_coeff) :: ctp1,ctp2
!\end{verbatim} %+
!  TYPE gtp_tpfun_as_coeff
! record for a TPFUN converted to coefficents without references to other funs
! Note ranges may increase when adding functions!!
!     double precision, dimension(:), allocatable :: tbreaks
!     double precision, dimension(:,:), allocatable :: coefs
!     integer, dimension(:,:), allocatable :: tpows
!  end type gtp_tpfun_as_coeff
! Max number of coefficients is maxnc, ranges is maxnr
! Previously I have used a maximum of 11 for any specific function ...   
! I am not sure I can just allocate bigger .... but I will try
   integer, parameter :: maxnc=15,maxnr=20
   double precision, parameter :: tenth=1.0D-1
   integer i1,i2,i3,k1,k2,nr3
   double precision tmax
   type(gtp_tpfun_as_coeff) :: ctp3
!
!   write(*,*)'In adjustranges, tmax: ',nr1,ctp1%tbreaks(nr1),ctp2%tbreaks(nr2)
   tmax=max(ctp1%tbreaks(nr1),ctp2%tbreaks(nr2))
   i1=1
   i2=1
   i3=0
   allocate(ctp3%tbreaks(maxnr))
   allocate(ctp3%coefs(maxnc,maxnr))
   allocate(ctp3%tpows(maxnc,maxnr))
   ctp3%tpows=-100
!----------------------------------------------------------------
!   write(*,79)'3Z adding and adjusting ranges: ',nr1,nr2
!79 format(a,2i2)
100 continue
   i3=i3+1
   if(i3.gt.maxnr) then
      write(*,*)'3Z too many ranges is summation function',i3
      gx%bmperr=4391; goto 1000
   endif
!   write(*,170)'3Z calling add1tp: ',i1,i2,i3,nr1,nr2
!170 format(a,10i5)
   call add1tpcoeffs(i1,ctp1,i2,ctp2,i3,ctp3)
   if(gx%bmperr.ne.0) goto 1000
   if(abs(ctp2%tbreaks(i2)-ctp1%tbreaks(i1)).lt.tenth) then
! both breakpoints the same!
!      write(*,180)'3Z same breakpoint',0,ctp1%tbreaks(i1),ctp2%tbreaks(i2)
!180   format(a,i3,2F10.2)
      ctp3%tbreaks(i3)=ctp2%tbreaks(i2)
      if(i2.lt.nr2) i2=i2+1
      if(i1.lt.nr1) i1=i1+1
   elseif(i1.eq.nr1 .and. i2.eq.nr2) then
!      write(*,180)'3Z breakpoint at max',0,tmax
      ctp3%tbreaks(i3)=tmax
   elseif(ctp2%tbreaks(i2).lt.ctp1%tbreaks(i1)) then
! we must create a breakpoint at the lowest tbreaks
      ctp3%tbreaks(i3)=ctp2%tbreaks(i2)
!     write(*,180)'3Z breakpoint in ctp2: ',i2,ctp2%tbreaks(i2),ctp1%tbreaks(i1)
      if(i2.lt.nr2) then
         i2=i2+1
      elseif(i1.lt.nr1) then
         i1=i1+1
      endif
   else
      ctp3%tbreaks(i3)=ctp1%tbreaks(i1)
!     write(*,180)'3Z breakpoint in ctp1: ',i1,ctp1%tbreaks(i1),ctp2%tbreaks(i2)
      if(i1.lt.nr1) then
         i1=i1+1
      elseif(i2.lt.nr2) then
         i2=i2+1
      endif
   endif
!   write(*,210)'3Z created ctp3 range: ',i3,ctp3%tbreaks(i3),tmax
!210 format(a,i3,2F9.2)
! How to know when we finished??
   if(abs(ctp3%tbreaks(i3)-tmax).gt.tenth) goto 100
!-------------------------------------------------------------
!   call tpwrite('jj',0,i3,ctp3)
   ctp1=ctp3
   nr1=i3
!   call tpwrite('hh',0,nr1,ctp1)
! I assume the arrays allocated for ctp3 will be deallocated automatically
1000 continue
   return
 end subroutine adjustranges

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine add1tpcoeffs
!\begin{verbatim} %-
 subroutine add1tpcoeffs(i1,ctp1,i2,ctp2,i3,ctp3)
! ctp3 is created with added coefficents from range i1 in ctp1 
! and range i2 in cp2 with same tpower.  Normally i3=1
   implicit none
   integer i1,i2,i3
   type(gtp_tpfun_as_coeff) :: ctp1,ctp2,ctp3
!\end{verbatim} %+
   integer, parameter :: maxnc=15
   integer j1,j2,j3,k1
! first copy ctp1 to ctp3.  Then add ctp3 coefficients with same powers
!   write(*,16)'3Z add1tp1: ',i1,i2,i3,ctp1%coefs(1,i1),ctp2%coefs(1,i2)
!16 format(a,3i3,2(1pe14.6))
!   write(*,17)(ctp1%tpows(j3,i1),j3=1,maxnc)
!   write(*,17)(ctp2%tpows(j3,i2),j3=1,maxnc)
!   write(*,17)(ctp3%tpows(j3,i3),j3=1,maxnc)
!17 format('3Z tpows: ',10i5)
   j3=0
!   call tpwrite('x0',0,i1,ctp1)
   do j1=1,maxnc
!      if(ctp1%tpows(j1,i1).gt.-100) then
         j3=j3+1
         ctp3%coefs(j3,i3)=ctp1%coefs(j1,i1)
         ctp3%tpows(j3,i3)=ctp1%tpows(j1,i1)
!      endif
   enddo
!   call tpwrite('x1',0,1,ctp3)
!   call tpwrite('x2',0,i2,ctp2)
!   write(*,*)'3Z no terms in ctp1?',j3
   f2: do j2=1,maxnc
      if(ctp2%tpows(j2,i2).gt.-100) then
         do j1=1,maxnc
            if(ctp2%tpows(j2,i2).eq.ctp3%tpows(j1,i3)) then
               ctp3%coefs(j1,i3)=ctp3%coefs(j1,i3)+ctp2%coefs(j2,i2)
               cycle f2
            endif
         enddo
! there is a t-power in ctp2 not already present in ctp3
!         write(*,*)'3Z have we managed to add?',j2,j3
         newpower: do j3=1,maxnc
            if(ctp3%tpows(j3,i3).le.-100) then
!               write(*,*)'3Z now we insert!',j3,ctp2%tpows(j2,i2),&
!                    ctp2%coefs(j2,i2)
               ctp3%coefs(j3,i3)=ctp2%coefs(j2,i2)
               ctp3%tpows(j3,i3)=ctp2%tpows(j2,i2)
               exit newpower 
            endif
         enddo newpower
      endif
   enddo f2
!   call tpwrite('x3',0,1,ctp3)
! that is all??
!   write(*,16)'3Z add1tp7: ',i1,i2,i3,ctp3%coefs(1,i3),ctp2%coefs(1,i2)
! the loops above may miss terms with same power ... suck
! check all terms in ctp3 
   do j1=1,maxnc
      if(ctp3%tpows(j1,i3).gt.-100) then
         do j2=j1+1,maxnc
            if(ctp3%tpows(j2,i3).eq.ctp3%tpows(j1,i3)) then
               ctp3%coefs(j1,i3)=ctp3%coefs(j1,i3)+ctp3%coefs(j2,i3)
               ctp3%tpows(j2,i3)=-100
               ctp3%coefs(j2,i3)=zero
            endif
         enddo
      endif
   enddo
!   call tpwrite('x4',0,1,ctp3)
!1000 continue
!   write(*,*)'3Z Exit add1tpcoefs'
   return
 end subroutine add1tpcoeffs

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine checkpowers
!\begin{verbatim} %-
 subroutine checkpowers(nc1,lfun,tpow1,npow,usedpow)
! check powers used in TP functions
! There can be several terms with same power ...
! nc1 is the maximal number of coefficients for each range (maxnc in fact)
   implicit none
   integer tpow1(*),nc1,lfun,usedpow(*),npow
!\end{verbatim} %+
! if these powers changes change also in sortceffs
!   integer, parameter :: fixpows(9)=[0,1,100,2,3,-1,7,-9,4]
   integer, parameter :: mmm=10
   integer, parameter :: fixpows(mmm)=[0,1,100,2,3,-1,7,-9,-2,-3]
! ANY CHANGE IN POWERS ALSO IN ... SORTCOEFFS
   integer i1,j1
   if(npow.eq.0) then
      do j1=1,nc1
         usedpow(j1)=-100
      enddo
!      npow=9
      npow=mmm
      do j1=1,npow
         usedpow(j1)=fixpows(j1)
      enddo
!      write(*,17)'3Z inititated usedpow: ',(usedpow(j1),j1=1,npow),lfun
!17    format(a,10i5,i3)
   endif
   loop1: do i1=1,nc1
      if(tpow1(i1).gt.-100) then
         do j1=1,npow
            if(tpow1(i1).eq.usedpow(j1)) cycle loop1
         enddo
! we have a non standard power
         npow=npow+1
         usedpow(npow)=tpow1(i1)
!         write(*,11)'3Z non-standard power: ',lfun,tpfuns(lfun)%symbol,&
!              npow,tpow1(i1)
!11       format(a,i5,2x,a,2i4)
      endif
   enddo loop1
1000 continue
   return
 end subroutine checkpowers

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine sortcoeffs
!\begin{verbatim} %-
 subroutine sortcoeffs(nc1,lfun,coeff1,tpow1)
! sort the coefficients in order power: 0 1 TlnT 2 3 -1; 7 -9 -2 other1 other2
!                                       1 2  3   4 5  6; 7  8  9   10    11
! other powers are 3, 4, -8 (meaning sqrt(t)) and maybe more
! tpowi is array giving the T power for coeffi,  101 means T*ln(T)
! ANY CHANGE OF POWERS MUST BE MADE ALSO IN ... CHECKPOWERS
! There can be several terms with same power ...
! nc1 is the maximal number of coefficients for each range (maxnc in fact)
   implicit none
   integer tpow1(*),nc1,lfun
   double precision coeff1(*)
!\end{verbatim}
   integer, parameter :: maxnc=15
   integer z1,i2,lastc,rare,free(3)
   double precision xxx,cord(maxnc)
   cord=zero
   lastc=0
   rare=0
   free=0
!   write(*,80)'3Z sort: ',tpfuns(lfun)%symbol,nc1,(tpow1(z1),z1=1,nc1)
   loop1: do z1=1,nc1
      if(tpow1(z1).eq.0) then
         cord(1)=cord(1)+coeff1(z1)
         if(lastc.lt.1) lastc=1
      elseif(tpow1(z1).eq.1) then
         cord(2)=cord(2)+coeff1(z1)
         if(lastc.lt.2) lastc=2
      elseif(tpow1(z1).eq.100) then
         cord(3)=cord(3)+coeff1(z1)
         if(lastc.lt.3) lastc=3
      elseif(tpow1(z1).eq.2) then
         cord(4)=cord(4)+coeff1(z1)
         if(lastc.lt.4) lastc=4
      elseif(tpow1(z1).eq.3) then
         cord(5)=cord(5)+coeff1(z1)
         if(lastc.lt.5) lastc=5
      elseif(tpow1(z1).eq.-1) then
         cord(6)=cord(6)+coeff1(z1)
         if(lastc.lt.6) lastc=6
      elseif(tpow1(z1).eq.7) then
! all powers from here are special ... if coeff(z1)=zero ignore on output
         cord(7)=cord(7)+coeff1(z1)
!         write(*,77)z1,7,tpfuns(lfun)%symbol
!77       format('3Z moving coefficient ',i2,' to ',i2,': ',a)
         coeff1(z1)=zero
         if(lastc.lt.7) lastc=7
      elseif(tpow1(z1).eq.-9) then
         cord(8)=cord(8)+coeff1(z1)
!         write(*,77)z1,8,tpfuns(lfun)%symbol
         coeff1(z1)=zero
         if(lastc.lt.8) lastc=8
      elseif(tpow1(z1).eq.-2) then
! it seems power -2 occors in the TAFID database
         cord(9)=cord(9)+coeff1(z1)
!         write(*,77)z1,9,tpfuns(lfun)%symbol
         coeff1(z1)=zero
         if(lastc.lt.9) lastc=9
      elseif(tpow1(z1).le.-100) then
! ignore this term
         continue
      elseif(coeff1(z1).ne.zero) then
! here tpow1(z1) cannot be -100: max 2 rare or unusual power like 3, 4, -8 ...
         if(free(1).eq.0 .or. tpow1(z1).eq.tpow1(10)) then
! store in unused or add to to same rare power position in position 10
!            write(*,90)tpfuns(lfun)%symbol,&
!                 lfun,10,z1,tpow1(z1),coeff1(z1)
!90          format('3Z function: ',a,' extra power: ',4i4,2x,1pe12.4)
            cord(10)=cord(10)+coeff1(z1)
            coeff1(z1)=zero
            tpow1(10)=tpow1(z1)
            free(1)=1
            if(lastc.lt.10) lastc=10
         elseif(free(2).eq.0 .or. tpow1(z1).eq.tpow1(11)) then
! store in unused or add to to same rare power position in position 11
! same special power in position 11
!            write(*,90)tpfuns(lfun)%symbol,&
!                 lfun,11,z1,tpow1(z1),coeff1(z1)
            cord(11)=cord(11)+coeff1(z1)
            coeff1(z1)=zero
            tpow1(11)=tpow1(z1)
            free(2)=1
            if(lastc.lt.11) lastc=11
         elseif(free(3).eq.0 .or. tpow1(z1).eq.tpow1(12)) then
! store in unused or add to to same rare power position in position 12
! same special power in position 12
!            write(*,90)tpfuns(lfun)%symbol,&
!                 lfun,11,z1,tpow1(z1),coeff1(z1)
            cord(12)=cord(12)+coeff1(z1)
            coeff1(z1)=zero
            tpow1(12)=tpow1(z1)
            free(3)=1
            if(lastc.lt.11) lastc=11
         else
! Too many rare powers in this expression
            write(*,89)tpfuns(lfun)%symbol,&
                 tpow1(10),tpow1(11),tpow1(12),tpow1(z1)
89          format('3Z Cannot handle four different rare powers: ',a,3i4)
            stop ' *** power problems!'
         endif
      endif
   enddo loop1
!   write(*,91)'3Z powers 1: ',(tpow1(z1),z1=1,lastc)
91 format(a,11i5)
! return coefficients in order
   do z1=1,9
      tpow1(z1)=-100
   enddo
   do z1=1,nc1
      coeff1(z1)=cord(z1)
   enddo
   tpow1(1)=0
   tpow1(2)=1
   tpow1(3)=100
   tpow1(4)=2
   tpow1(5)=3
   tpow1(6)=-1
   tpow1(7)=7
   tpow1(8)=-9
   tpow1(9)=-2
! latsc is the last used power position   
   if(lastc.lt.10) tpow1(10)=-100
   if(lastc.lt.11) tpow1(11)=-100
   if(lastc.lt.12) tpow1(12)=-100
   tpow1(13)=-100; tpow1(13)=-100; tpow1(13)=-100
!   tpow1(10), 11 and 12 are free
!   write(*,91)'3Z powers 2: ',(tpow1(z1),z1=1,lastc)
!   write(*,80)'3Z sorted: ',tpfuns(lfun)%symbol,nc1,(tpow1(z1),z1=1,lastc)
!80 format(a,1x,a,': ',i3,11i5)
! tpow1(11) keep its value.  No provision for more than one extra power!!
!1000 continue
   return
 end subroutine sortcoeffs

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

! END MODULE TPFUNLIB

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

