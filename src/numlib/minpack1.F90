! This extract from MINPACK contains:
! LMDIF1: least square routine
! HYBRD1: solving systems of non-linear equations
! and some support routines
! and calfun calling OC for assessment_calfun for optimization
!
!MODULE LIBOCEQPLUS
!
MODULE MINPACK
!
!  use liboceq
!  use minpack2
!
!
!  Minpack Copyright Notice (1999) University of Chicago.  All rights reserved
! 
! Redistribution and use in source and binary forms, with or
! without modification, are permitted provided that the
! following conditions are met:
! 
! 1. Redistributions of source code must retain the above
! copyright notice, this list of conditions and the following
! disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above
! copyright notice, this list of conditions and the following
! disclaimer in the documentation and/or other materials
! provided with the distribution.
! 
! 3. The end-user documentation included with the
! redistribution, if any, must include the following
! acknowledgment:
! 
!    "This product includes software developed by the
!    University of Chicago, as Operator of Argonne National
!    Laboratory.
! 
! Alternately, this acknowledgment may appear in the software
! itself, if and wherever such third-party acknowledgments
! normally appear.
! 
! 4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
! WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
! UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
! THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
! OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
! OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
! USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
! THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
! DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
! UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
! BE CORRECTED.
! 
! 5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
! HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
! ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
! INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
! ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
! PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
! SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
! (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
! EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
! POSSIBILITY OF SUCH LOSS OR DAMAGES.
!
  implicit none
  double precision, parameter, private :: zero=0.0d0
!
contains
!
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  subroutine lmdif1(fcn,m,n,x,fvec,tol,info,nfev,iwa,wa,lwa,fjac,err0)
! call modified by Bo Sundman 2017, modified again 2018 to include fcn
! nfev is number of calls to fcn
!  subroutine lmdif1(    m,n,x,fvec,tol,info,nfev,iwa,wa,lwa,fjac,err0)
! original:
!  subroutine lmdif1(fcn,m,n,x,fvec,tol,info,     iwa,wa,lwa)
    implicit none
    integer m,n,info,lwa
    integer iwa(n)
    double precision tol
    double precision x(n),fvec(m),wa(lwa),fjac(m,*)
    external fcn
!     **********
!
!     subroutine lmdif1
!
!     the purpose of lmdif1 is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of the
!     levenberg-marquardt algorithm. this is done by using the more
!     general least-squares solver lmdif. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!     the subroutine statement is
!
!       subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
! modified to include iterations: subroutine fcn(m,n,x,fvec,iflag,nfev)
!         integer m,n,iflag,nfev
!         double precision x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmdif1.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates either that the relative
!         error in the sum of squares is at most tol or that
!         the relative error between x and the solution is at
!         most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  algorithm estimates that the relative error
!                   in the sum of squares is at most tol.
!
!         info = 2  algorithm estimates that the relative error
!                   between x and the solution is at most tol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  fvec is orthogonal to the columns of the
!                   jacobian to machine precision.
!
!         info = 5  number of calls to fcn has reached or
!                   exceeded 200*(n+1).
!
!         info = 6  tol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  tol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!       iwa is an integer work array of length n.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than
!         m*n+5*n+m.
!
!       fjac added to calculate relative standard deviation (SD)
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... lmdif
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
    integer maxfev,mode,mp5n,nfev,nprint
!      double precision epsfcn,factor,ftol,gtol,xtol,zero
! err0 contains intitial sum of error and last sum of errors
    double precision epsfcn,ftol,gtol,xtol,err0(*)
!      data factor,zero /1.0d2,0.0d0/
! zero already defined globally
    double precision :: factor=1.0D2,zero=0.0D0
! number of iterations passed through infor
    maxfev=info
    info = 0
!    write(*,*)'in lmdif1',maxfev,m,n
!
!     check the input parameters for errors.
!
    if (n .le. 0 .or. m .lt. n .or. tol .lt. zero &
         .or. lwa .lt. m*n + 5*n + m) then
       write(*,*)'Illegal call of lmdif1'
       go to 10
    endif
!
!     call lmdif.
!
    info=0
! several of these moved to lmdif ... as well as allocating workspace
    ftol = tol
    xtol = tol
    gtol = zero
    epsfcn = zero
    mode = 1
! This controls output during optimization, output must be added to calfun
!    nprint = 0
    nprint=1
    mp5n = m + 5*n
!    write(*,*)'Calling lmdif1',maxfev
!    call lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,wa(1), &
!                mode,factor,nprint,info,nfev,wa(mp5n+1),m,iwa, &
!                wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
! remove fcn and reduce number of arguments and linker chokes ...
!    write(*,*)'minpack: lmdif1 call lmdif',m,n
    call lmdif(fcn,m,n,x,fvec,tol,maxfev, &
                mode,factor,nprint,info,nfev,fjac,iwa,err0)
!    call lmdif(m,n,x,fvec,tol,maxfev, &
!                mode,factor,nprint,info,nfev,fjac,iwa,err0)
    if (info .eq. 8) info = 4
!    write(*,*)'Return from lmdif with info= ',info
10  continue
    return
!
!     last card of subroutine lmdif1.
!
  end subroutine lmdif1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!  subroutine fdjac2(m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
  subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
    implicit none
    integer m,n,ldfjac,iflag
    double precision epsfcn
    double precision x(n),fvec(m),fjac(ldfjac,n),wa(m)
!     **********
!
!     subroutine fdjac2
!
!     this subroutine computes a forward-difference approximation
!     to the m by n jacobian matrix associated with a specified
!     problem of m functions in n variables.
!
!     the subroutine statement is
!
!       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
! modified to include iterations: subroutine fcn(m,n,x,fvec,iflag,niter)
!         integer m,n,iflag,niter
!         double precision x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of fdjac2.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an input array of length n.
!
!       fvec is an input array of length m which must contain the
!         functions evaluated at x.
!
!       fjac is an output m by n array which contains the
!         approximation to the jacobian matrix evaluated at x.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       iflag is an integer variable which can be used to terminate
!         the execution of fdjac2. see description of fcn.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       wa is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar
!
!       fortran-supplied ... dabs,dmax1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
    integer i,j
    double precision eps,epsmch,h,temp,zero
!    integer, parameter :: niter=0
! added as argumnet added to fcn and it may be incremented in fcn ...
    integer :: niter=0
! missing external declaration
    external fcn
!    double precision dpmpar
!    data zero /0.0d0/
!
!     epsmch is the machine precision.
!
    epsmch = dpmpar(1)
!
    eps = dsqrt(dmax1(epsfcn,epsmch))
!    write(*,17)'epsis: ',0,epsfcn,epsmch,eps
!    eps=1.0D-4
!      do 20 j = 1, n
    do j = 1, n
       temp = x(j)
       h = eps*dabs(temp)
       if (h .eq. zero) h = eps
!       write(*,17)'In fdjac2: ',j,temp,h,eps
!17     format(a,i2,6(1pe12.4))
       x(j) = temp + h
!       call fcn(m,n,x,wa,iflag)    <<<<<<<<<<<< original
       call fcn(m,n,x,wa,iflag,niter)
!       call calfun(m,n,x,wa,iflag,niter)
       if (iflag .lt. 0) go to 30
       x(j) = temp
!         do 10 i = 1, m
       do i = 1, m
          fjac(i,j) = (wa(i) - fvec(i))/h
       enddo
    enddo
!   10       continue
!   20    continue
30  continue
    return
!
!     last card of subroutine fdjac2.
!
  end subroutine fdjac2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!  subroutine lmdif(m,n,x,fvec,xtol,maxfev, &
  subroutine lmdif(fcn,m,n,x,fvec,xtol,maxfev, &
       mode,factor,nprint,info,nfev,fjac,ipvt,err0)
! removed arguments as linker chokes ...
!  subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,diag, &
!       mode,factor,nprint,info,nfev,fjac,ldfjac, &
!       ipvt,qtf,wa1,wa2,wa3,wa4)
    implicit none
    integer m,n,maxfev,mode,nprint,info,nfev,ldfjac
    integer ipvt(n)
    double precision ftol,xtol,gtol,epsfcn,factor,err0(*)
    double precision x(n),fvec(m),fjac(m,*)
!    double precision x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n), &
!    double precision x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n), &
!         wa1(n),wa2(n),wa3(n),wa4(m)
    external fcn
!     **********
    double precision, dimension(:), allocatable :: diag,qtf,wa1,wa2,wa3,wa4
!    double precision, dimension(:,:), allocatable :: fjac
!
!     subroutine lmdif
!
!     the purpose of lmdif is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!     the subroutine statement is
!
!       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
!                        diag,mode,factor,nprint,info,nfev,fjac,
!                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
! modified to iiterations: subroutine fcn(m,n,x,fvec,iflag,niter)
!         integer m,n,iflag,niter
!         double precision x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmdif.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn is at least
!         maxfev by the end of an iteration.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  both actual and predicted relative reductions
!                   in the sum of squares are at most ftol.
!
!         info = 2  relative error between two consecutive iterates
!                   is at most xtol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  the cosine of the angle between fvec and any
!                   column of the jacobian is at most gtol in
!                   absolute value.
!
!         info = 5  number of calls to fcn has reached or
!                   exceeded maxfev.
!
!         info = 6  ftol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  xtol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!         info = 8  gtol is too small. fvec is orthogonal to the
!                   columns of the jacobian to machine precision.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
!       wa1, wa2, and wa3 are work arrays of length n.
!
!       wa4 is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
    integer i,iflag,iter,j,l
    double precision actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm, &
         par,pnorm,prered,ratio, &
         sum,temp,temp1,temp2,xnorm,bosum
! removed variables one line above 
!         one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio, &
! if the functions dpmpar and enorm are declared here strange link error ...
!    double precision dpmpar,enorm
!    data one,p1,p5,p25,p75,p0001,zero &
!         /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
    double precision :: one=1.0D0
    double precision :: p1=1.0D-1,p5=5.0D-1,p25=2.5D-1,p75=7.5D-1,p0001=1.0D-4
!
!    write(*,*)'minpack: in lmdif A: ',n,m
! replace removed arguments
    ldfjac=m
    ftol=xtol
    gtol=zero
    epsfcn=zero
    allocate(diag(n)) 
    allocate(qtf(n)) 
    allocate(wa1(n))
    allocate(wa2(n))
    allocate(wa3(n))
    allocate(wa4(m))
! now included in call
!    allocate(fjac(ldfjac,n))
!
!     epsmch is the machine precision.
!
    epsmch = dpmpar(1)
!
    info = 0
    iflag = 0
    nfev = 0
!
!     check the input parameters for errors.
!
! modified to run once if maxfev=0 (dry run)
!    write(*,*)'In lmdif C: maxfev=',maxfev,m,n
    if (n .le. 0 .or. m .lt. n .or. ldfjac .lt. m &
         .or. ftol .lt. zero .or. xtol .lt. zero .or. gtol .lt. zero &
         .or. maxfev .lt. 0 .or. factor .le. zero) go to 300
    if (mode .ne. 2) go to 20
    do j = 1, n
       if (diag(j) .le. zero) go to 300
    enddo
!   10    continue
20  continue
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
    iflag = 1
!    call fcn(m,n,x,fvec,iflag)             <<<<<<<<<<< original
!    call calfun(m,n,x,fvec,iflag,nfev)
!    write(*,*)'minpack: lmdif call fcn 1: ',n,m
    call fcn(m,n,x,fvec,iflag,nfev)
! calculate intial sum of errors
!    write(*,*)'lmdif back from calfun',nfev
    bosum=zero
    do j=1,m
       bosum=bosum+fvec(j)**2
    enddo
    err0(1)=bosum
!--------------------------
    nfev = 1
    if (iflag .lt. 0) go to 300
    if(maxfev .eq. 0) goto 300
    fnorm = enorm(m,fvec)
!
!     initialize levenberg-marquardt parameter and iteration counter.
!
    par = zero
    iter = 1
!
!     beginning of the outer loop.
!
30  continue
!
!        calculate the jacobian matrix.
!
    iflag = 2
!    write(*,*)'minpack: lmdif call fdjac2 1: ',n,m
    call fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa4)
!    call fdjac2(m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa4)
    nfev = nfev + n
    if (iflag .lt. 0) go to 300
!
!        if requested, call fcn to enable printing of iterates.
!
    if (nprint .le. 0) go to 40
    iflag = 0
    if (mod(iter-1,nprint) .eq. 0) then
!       call calfun(m,n,x,fvec,iflag,nfev)
!       write(*,*)'minpack: lmdif call fcn 3: ',n,m
       call fcn(m,n,x,fvec,iflag,nfev)
    endif
!    if (mod(iter-1,nprint) .eq. 0) call fcn(m,n,x,fvec,iflag)
    if (iflag .lt. 0) go to 300
40  continue
!
!        compute the qr factorization of the jacobian.
!
    call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
    if (iter .ne. 1) go to 80
    if (mode .eq. 2) go to 60
!    do 50 j = 1, n
    do j = 1, n
       diag(j) = wa2(j)
       if (wa2(j) .eq. zero) diag(j) = one
    enddo
!   50       continue
60  continue
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
    do j = 1, n
       wa3(j) = diag(j)*x(j)
    enddo
!   70       continue
    xnorm = enorm(n,wa3)
    delta = factor*xnorm
    if (delta .eq. zero) delta = factor
80  continue
!
!        form (q transpose)*fvec and store the first n components in
!        qtf.
!
    do i = 1, m
       wa4(i) = fvec(i)
    enddo
!90  continue
!    do 130 j = 1, n
    do j = 1, n
       if (fjac(j,j) .eq. zero) go to 120
       sum = zero
!       do 100 i = j, m
       do i = j, m
          sum = sum + fjac(i,j)*wa4(i)
       enddo
!100       continue
       temp = -sum/fjac(j,j)
       do i = j, m
          wa4(i) = wa4(i) + fjac(i,j)*temp
       enddo
!  110          continue
120    continue
       fjac(j,j) = wa1(j)
       qtf(j) = wa4(j)
    enddo
!130    continue
!
!        compute the norm of the scaled gradient.
!
    gnorm = zero
    if (fnorm .eq. zero) go to 170
!    do 160 j = 1, n
    do j = 1, n
       l = ipvt(j)
       if (wa2(l) .eq. zero) go to 150
       sum = zero
       do i = 1, j
          sum = sum + fjac(i,j)*(qtf(i)/fnorm)
       enddo
!  140          continue
       gnorm = dmax1(gnorm,dabs(sum/wa2(l)))
150    continue
    enddo
!160    continue
170    continue
!
!        test for convergence of the gradient norm.
!
    if (gnorm .le. gtol) info = 4
    if (info .ne. 0) go to 300
!
!        rescale if necessary.
!
    if (mode .eq. 2) go to 190
    do j = 1, n
       diag(j) = dmax1(diag(j),wa2(j))
    enddo
!180         continue
190 continue
!
!        beginning of the inner loop.
!
200 continue
!
!           determine the levenberg-marquardt parameter.
!
    call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2, &
         wa3,wa4)
!
!           store the direction p and x + p. calculate the norm of p.
!
!    do 210 j = 1, n
    do j = 1, n
       wa1(j) = -wa1(j)
       wa2(j) = x(j) + wa1(j)
       wa3(j) = diag(j)*wa1(j)
    enddo
!210 continue
    pnorm = enorm(n,wa3)
!
!           on the first iteration, adjust the initial step bound.
!
    if (iter .eq. 1) delta = dmin1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
    iflag = 1
!    call fcn(m,n,wa2,wa4,iflag)
!    call calfun(m,n,wa2,wa4,iflag,nfev)
!    write(*,*)'minpack: lmdif call fcn 3: ',n,m
    call fcn(m,n,wa2,wa4,iflag,nfev)
    nfev = nfev + 1
    if (iflag .lt. 0) go to 300
    fnorm1 = enorm(m,wa4)
!
!           compute the scaled actual reduction.
!
    actred = -one
    if (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction and
!           the scaled directional derivative.
!
!    do 230 j = 1, n
    do j = 1, n
       wa3(j) = zero
       l = ipvt(j)
       temp = wa1(l)
       do i = 1, j
          wa3(i) = wa3(i) + fjac(i,j)*temp
       enddo
    enddo
!220       continue
!230    continue
    temp1 = enorm(n,wa3)/fnorm
    temp2 = (dsqrt(par)*pnorm)/fnorm
    prered = temp1**2 + temp2**2/p5
    dirder = -(temp1**2 + temp2**2)
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
    ratio = zero
    if (prered .ne. zero) ratio = actred/prered
!
!           update the step bound.
!
    if (ratio .gt. p25) go to 240
    if (actred .ge. zero) temp = p5
    if (actred .lt. zero) &
         temp = p5*dirder/(dirder + p5*actred)
    if (p1*fnorm1 .ge. fnorm .or. temp .lt. p1) temp = p1
    delta = temp*dmin1(delta,pnorm/p1)
    par = par/temp
    go to 260
240 continue
    if (par .ne. zero .and. ratio .lt. p75) go to 250
    delta = pnorm/p5
    par = p5*par
250 continue
260 continue
!
!           test for successful iteration.
!
    if (ratio .lt. p0001) go to 290
!
!           successful iteration. update x, fvec, and their norms.
!
    do j = 1, n
       x(j) = wa2(j)
       wa2(j) = diag(j)*x(j)
    enddo
!270    continue
    do i = 1, m
       fvec(i) = wa4(i)
    enddo
!280    continue
    xnorm = enorm(n,wa2)
    fnorm = fnorm1
    iter = iter + 1
290 continue
!
!           tests for convergence.
!
    if (dabs(actred) .le. ftol .and. prered .le. ftol &
         .and. p5*ratio .le. one) info = 1
    if (delta .le. xtol*xnorm) info = 2
    if (dabs(actred) .le. ftol .and. prered .le. ftol &
         .and. p5*ratio .le. one .and. info .eq. 2) info = 3
    if (info .ne. 0) go to 300
!
!           tests for termination and stringent tolerances.
!
    if (nfev .ge. maxfev) info = 5
    if (dabs(actred) .le. epsmch .and. prered .le. epsmch &
         .and. p5*ratio .le. one) info = 6
    if (delta .le. epsmch*xnorm) info = 7
    if (gnorm .le. epsmch) info = 8
    if (info .ne. 0) go to 300
!
!           end of the inner loop. repeat if iteration unsuccessful.
!
    if (ratio .lt. p0001) go to 200
!
!        end of the outer loop.
!
    go to 30
300 continue
!
!     termination, either normal or user imposed.
!
    if (iflag .lt. 0) info = iflag
    iflag = 0
!    write(*,*)'minpack: lmdif call fcn 4: ',n,m,info,maxfev
    if(maxfev.gt.0) then
       if (nprint .gt. 0) call fcn(m,n,x,fvec,iflag)
!    if (nprint .gt. 0) call calfun(m,n,x,fvec,iflag,-nfev)
    else
! Add that calfun called once if maxfev=0 to calculate all errors
!    if (maxfev .eq. 0) call calfun(m,n,x,fvec,1,0)
       call fcn(m,n,x,fvec,1,0)
    endif
!    write(*,*)'minpack: lmdif call fcn 5: ',n,m,maxfev
    return
!
!     last card of subroutine lmdif.
!
  end subroutine lmdif

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1, &
       wa2)
    implicit none
    integer n,ldr
    integer ipvt(n)
    double precision delta,par
    double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa1(n), &
         wa2(n)
!     **********
!
!     subroutine lmpar
!
!     given an m by n matrix a, an n by n nonsingular diagonal
!     matrix d, an m-vector b, and a positive number delta,
!     the problem is to determine a value for the parameter
!     par such that if x solves the system
!
!           a*x = b ,     sqrt(par)*d*x = 0 ,
!
!     in the least squares sense, and dxnorm is the euclidean
!     norm of d*x, then either par is zero and
!
!           (dxnorm-delta) .le. 0.1*delta ,
!
!     or par is positive and
!
!           abs(dxnorm-delta) .le. 0.1*delta .
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization, with column pivoting, of a. that is, if
!     a*p = q*r, where p is a permutation matrix, q has orthogonal
!     columns, and r is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then lmpar expects
!     the full upper triangle of r, the permutation matrix p,
!     and the first n components of (q transpose)*b. on output
!     lmpar also provides an upper triangular matrix s such that
!
!            t   t                   t
!           p *(a *a + par*d*d)*p = s *s .
!
!     s is employed within lmpar and may be of separate interest.
!
!     only a few iterations are generally needed for convergence
!     of the algorithm. if, however, the limit of 10 iterations
!     is reached, then the output par will contain the best
!     value obtained so far.
!
!     the subroutine statement is
!
!       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
!                        wa1,wa2)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. on input the full upper triangle
!         must contain the full upper triangle of the matrix r.
!         on output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix s.
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       ipvt is an integer input array of length n which defines the
!         permutation matrix p such that a*p = q*r. column j of p
!         is column ipvt(j) of the identity matrix.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.
!
!       delta is a positive input variable which specifies an upper
!         bound on the euclidean norm of d*x.
!
!       par is a nonnegative variable. on input par contains an
!         initial estimate of the levenberg-marquardt parameter.
!         on output par contains the final estimate.
!
!       x is an output array of length n which contains the least
!         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
!         for the output par.
!
!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix s.
!
!       wa1 and wa2 are work arrays of length n.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar,enorm,qrsolv
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
    integer i,iter,j,jm1,jp1,k,l,nsing
    double precision dxnorm,dwarf,fp,gnorm,parc,parl,paru,p1,p001, &
         sum,temp,zero
!    double precision dpmpar,enorm
    data p1,p001,zero /1.0d-1,1.0d-3,0.0d0/
!
!     dwarf is the smallest positive magnitude.
!
    dwarf = dpmpar(2)
!
!     compute and store in x the gauss-newton direction. if the
!     jacobian is rank-deficient, obtain a least squares solution.
!
    nsing = n
    do j = 1, n
       wa1(j) = qtb(j)
       if (r(j,j) .eq. zero .and. nsing .eq. n) nsing = j - 1
       if (nsing .lt. n) wa1(j) = zero
    enddo
!   10    continue
    if (nsing .lt. 1) go to 50
!      do 40 k = 1, nsing
    do k = 1, nsing
       j = nsing - k + 1
       wa1(j) = wa1(j)/r(j,j)
       temp = wa1(j)
       jm1 = j - 1
       if (jm1 .lt. 1) go to 30
       do i = 1, jm1
          wa1(i) = wa1(i) - r(i,j)*temp
       enddo
!20          continue
30     continue
    enddo
!40     continue
50  continue
    do j = 1, n
       l = ipvt(j)
       x(l) = wa1(j)
    enddo
!60     continue
!
!     initialize the iteration counter.
!     evaluate the function at the origin, and test
!     for acceptance of the gauss-newton direction.
!
    iter = 0
    do j = 1, n
       wa2(j) = diag(j)*x(j)
    enddo
!70     continue
    dxnorm = enorm(n,wa2)
    fp = dxnorm - delta
    if (fp .le. p1*delta) go to 220
!
!     if the jacobian is not rank deficient, the newton
!     step provides a lower bound, parl, for the zero of
!     the function. otherwise set this bound to zero.
!
    parl = zero
    if (nsing .lt. n) go to 120
    do j = 1, n
       l = ipvt(j)
       wa1(j) = diag(l)*(wa2(l)/dxnorm)
    enddo
!80     continue
!    do 110 j = 1, n
    do j = 1, n
       sum = zero
       jm1 = j - 1
       if (jm1 .lt. 1) go to 100
       do i = 1, jm1
          sum = sum + r(i,j)*wa1(i)
       enddo
!90        continue
100    continue
       wa1(j) = (wa1(j) - sum)/r(j,j)
    enddo
!110 continue
    temp = enorm(n,wa1)
    parl = ((fp/delta)/temp)/temp
120 continue
!
!     calculate an upper bound, paru, for the zero of the function.
!
!      do 140 j = 1, n
    do j = 1, n
       sum = zero
       do i = 1, j
          sum = sum + r(i,j)*qtb(i)
       enddo
!130    continue
       l = ipvt(j)
       wa1(j) = sum/diag(l)
    enddo
!140 continue
    gnorm = enorm(n,wa1)
    paru = gnorm/delta
    if (paru .eq. zero) paru = dwarf/dmin1(delta,p1)
!
!     if the input par lies outside of the interval (parl,paru),
!     set par to the closer endpoint.
!
    par = dmax1(par,parl)
    par = dmin1(par,paru)
    if (par .eq. zero) par = gnorm/dxnorm
!
!     beginning of an iteration.
!
150 continue
    iter = iter + 1
!
!        evaluate the function at the current value of par.
!
    if (par .eq. zero) par = dmax1(dwarf,p001*paru)
    temp = dsqrt(par)
    do j = 1, n
       wa1(j) = temp*diag(j)
    enddo
!160 continue
    call qrsolv(n,r,ldr,ipvt,wa1,qtb,x,sdiag,wa2)
    do j = 1, n
       wa2(j) = diag(j)*x(j)
    enddo
!170 continue
    dxnorm = enorm(n,wa2)
    temp = fp
    fp = dxnorm - delta
!
!        if the function is small enough, accept the current value
!        of par. also test for the exceptional cases where parl
!        is zero or the number of iterations has reached 10.
!
    if (dabs(fp) .le. p1*delta &
         .or. parl .eq. zero .and. fp .le. temp &
         .and. temp .lt. zero .or. iter .eq. 10) go to 220
!
!        compute the newton correction.
!
    do j = 1, n
       l = ipvt(j)
       wa1(j) = diag(l)*(wa2(l)/dxnorm)
    enddo
!180 continue
!         do 210 j = 1, n
    do j = 1, n
       wa1(j) = wa1(j)/sdiag(j)
       temp = wa1(j)
       jp1 = j + 1
       if (n .lt. jp1) go to 200
       do i = jp1, n
          wa1(i) = wa1(i) - r(i,j)*temp
       enddo
!190    continue
200    continue
    enddo
!210 continue
    temp = enorm(n,wa1)
    parc = ((fp/delta)/temp)/temp
!
!        depending on the sign of the function, update parl or paru.
!
    if (fp .gt. zero) parl = dmax1(parl,par)
    if (fp .lt. zero) paru = dmin1(paru,par)
!
!        compute an improved estimate for par.
!
    par = dmax1(parl,par+parc)
!
!        end of an iteration.
!
    go to 150
220 continue
!
!     termination.
!
    if (iter .eq. 0) par = zero
    return
!
!     last card of subroutine lmpar.
!
  end subroutine lmpar

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
    implicit none
    integer n,ldr
    integer ipvt(n)
    double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n)
!     **********
!
!     subroutine qrsolv
!
!     given an m by n matrix a, an n by n diagonal matrix d,
!     and an m-vector b, the problem is to determine an x which
!     solves the system
!
!           a*x = b ,     d*x = 0 ,
!
!     in the least squares sense.
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization, with column pivoting, of a. that is, if
!     a*p = q*r, where p is a permutation matrix, q has orthogonal
!     columns, and r is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then qrsolv expects
!     the full upper triangle of r, the permutation matrix p,
!     and the first n components of (q transpose)*b. the system
!     a*x = b, d*x = 0, is then equivalent to
!
!                  t       t
!           r*z = q *b ,  p *d*p*z = 0 ,
!
!     where x = p*z. if this system does not have full rank,
!     then a least squares solution is obtained. on output qrsolv
!     also provides an upper triangular matrix s such that
!
!            t   t               t
!           p *(a *a + d*d)*p = s *s .
!
!     s is computed within qrsolv and may be of separate interest.
!
!     the subroutine statement is
!
!       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. on input the full upper triangle
!         must contain the full upper triangle of the matrix r.
!         on output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix s.
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       ipvt is an integer input array of length n which defines the
!         permutation matrix p such that a*p = q*r. column j of p
!         is column ipvt(j) of the identity matrix.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.
!
!       x is an output array of length n which contains the least
!         squares solution of the system a*x = b, d*x = 0.
!
!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix s.
!
!       wa is a work array of length n.
!
!     subprograms called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
    integer i,j,jp1,k,kp1,l,nsing
    double precision cos,cotan,p5,p25,qtbpj,sin,sum,tan,temp,zero
    data p5,p25,zero /5.0d-1,2.5d-1,0.0d0/
!
!     copy r and (q transpose)*b to preserve input and initialize s.
!     in particular, save the diagonal elements of r in x.
!
!    do 20 j = 1, n
    do j = 1, n
       do i = j, n
          r(i,j) = r(j,i)
       enddo
!10     continue
       x(j) = r(j,j)
       wa(j) = qtb(j)
    enddo
!20  continue
!
!     eliminate the diagonal matrix d using a givens rotation.
!
!      do 100 j = 1, n
    do j = 1, n
!
!        prepare the row of d to be eliminated, locating the
!        diagonal element using p from the qr factorization.
!
       l = ipvt(j)
       if (diag(l) .eq. zero) go to 90
       do k = j, n
          sdiag(k) = zero
       enddo
!30     continue
       sdiag(j) = diag(l)
!
!        the transformations to eliminate the row of d
!        modify only a single element of (q transpose)*b
!        beyond the first n, which is initially zero.
!
       qtbpj = zero
!       do 80 k = j, n
       do k = j, n
!
!           determine a givens rotation which eliminates the
!           appropriate element in the current row of d.
!
          if (sdiag(k) .eq. zero) go to 70
          if (dabs(r(k,k)) .ge. dabs(sdiag(k))) go to 40
          cotan = r(k,k)/sdiag(k)
          sin = p5/dsqrt(p25+p25*cotan**2)
          cos = sin*cotan
          go to 50
40        continue
          tan = sdiag(k)/r(k,k)
          cos = p5/dsqrt(p25+p25*tan**2)
          sin = cos*tan
50        continue
!
!           compute the modified diagonal element of r and
!           the modified element of ((q transpose)*b,0).
!
          r(k,k) = cos*r(k,k) + sin*sdiag(k)
          temp = cos*wa(k) + sin*qtbpj
          qtbpj = -sin*wa(k) + cos*qtbpj
          wa(k) = temp
!
!           accumulate the tranformation in the row of s.
!
          kp1 = k + 1
          if (n .lt. kp1) go to 70
          do i = kp1, n
             temp = cos*r(i,k) + sin*sdiag(i)
             sdiag(i) = -sin*r(i,k) + cos*sdiag(i)
             r(i,k) = temp
          enddo
!60        continue
70        continue
       enddo
!80     continue
90     continue
!
!        store the diagonal element of s and restore
!        the corresponding diagonal element of r.
!
       sdiag(j) = r(j,j)
       r(j,j) = x(j)
    enddo
100 continue
!
!     solve the triangular system for z. if the system is
!     singular, then obtain a least squares solution.
!
    nsing = n
    do j = 1, n
       if (sdiag(j) .eq. zero .and. nsing .eq. n) nsing = j - 1
       if (nsing .lt. n) wa(j) = zero
    enddo
!110 continue
    if (nsing .lt. 1) go to 150
!    do 140 k = 1, nsing
    do k = 1, nsing
       j = nsing - k + 1
       sum = zero
       jp1 = j + 1
       if (nsing .lt. jp1) go to 130
       do i = jp1, nsing
          sum = sum + r(i,j)*wa(i)
       enddo
!120    continue
130    continue
       wa(j) = (wa(j) - sum)/sdiag(j)
    enddo
!140 continue
150 continue
!
!     permute the components of z back to components of x.
!
    do j = 1, n
       l = ipvt(j)
       x(l) = wa(j)
    enddo
!160 continue
    return
!
!     last card of subroutine qrsolv.
!
  end subroutine qrsolv

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\


  subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
    implicit none
    integer n,info,lwa
    double precision tol
    double precision x(n),fvec(n),wa(lwa)
    external fcn
!     **********
!
!     subroutine hybrd1
!
!     the purpose of hybrd1 is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. this is done by using the
!     more general nonlinear equation solver hybrd. the user
!     must provide a subroutine which calculates the functions.
!     the jacobian is then calculated by a forward-difference
!     approximation.
!
!     the subroutine statement is
!
!       subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,iflag)
!         integer n,iflag
!         double precision x(n),fvec(n)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrd1.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length n which contains
!         the functions evaluated at the output x.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates that the relative error
!         between x and the solution is at most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0   improper input parameters.
!
!         info = 1   algorithm estimates that the relative error
!                    between x and the solution is at most tol.
!
!         info = 2   number of calls to fcn has reached or exceeded
!                    200*(n+1).
!
!         info = 3   tol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than
!         (n*(3*n+13))/2.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... hybrd
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
    integer index,j,lr,maxfev,ml,mode,mu,nfev,nprint
    double precision epsfcn,factor,one,xtol,zero
    data factor,one,zero /1.0d2,1.0d0,0.0d0/
    info = 0
!
!     check the input parameters for errors.
!
    if (n .le. 0 .or. tol .lt. zero .or. lwa .lt. (n*(3*n + 13))/2) &
         go to 20
!
!     call hybrd.
!
    maxfev = 200*(n + 1)
    xtol = tol
    ml = n - 1
    mu = n - 1
    epsfcn = zero
    mode = 2
    do j = 1, n
       wa(j) = one
    enddo
    nprint = 0
    lr = (n*(n + 1))/2
    index = 6*n + lr
    call hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,wa(1),mode, &
         factor,nprint,info,nfev,wa(index+1),n,wa(6*n+1),lr, &
         wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
    if (info .eq. 5) info = 4
20  continue
    return
!
!     last card of subroutine hybrd1.
!                       
  end subroutine hybrd1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag, &
       mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr, &
       qtf,wa1,wa2,wa3,wa4)
    implicit none
    integer n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr
    double precision xtol,epsfcn,factor
    double precision x(n),fvec(n),diag(n),fjac(ldfjac,n),r(lr), &
         qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)
    external fcn
!     **********
!
!     subroutine hybrd
!
!     the purpose of hybrd is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!     the subroutine statement is
!
!       subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,
!                        diag,mode,factor,nprint,info,nfev,fjac,
!                        ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,iflag)
!         integer n,iflag
!         double precision x(n),fvec(n)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrd.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length n which contains
!         the functions evaluated at the output x.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn is at least maxfev
!         by the end of an iteration.
!
!       ml is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         ml to at least n - 1.
!
!       mu is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         mu to at least n - 1.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0   improper input parameters.
!
!         info = 1   relative error between two consecutive iterates
!                    is at most xtol.
!
!         info = 2   number of calls to fcn has reached or exceeded
!                    maxfev.
!
!         info = 3   xtol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    five jacobian evaluations.
!
!         info = 5   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    ten iterations.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn.
!
!       fjac is an output n by n array which contains the
!         orthogonal matrix q produced by the qr factorization
!         of the final approximate jacobian.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       r is an output array of length lr which contains the
!         upper triangular matrix produced by the qr factorization
!         of the final approximate jacobian, stored rowwise.
!
!       lr is a positive integer input variable not less than
!         (n*(n+1))/2.
!
!       qtf is an output array of length n which contains
!         the vector (q transpose)*fvec.
!
!       wa1, wa2, wa3, and wa4 are work arrays of length n.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dogleg,dpmpar,enorm,fdjac1,
!                            qform,qrfac,r1mpyq,r1updt
!
!       fortran-supplied ... dabs,dmax1,dmin1,min0,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
    integer i,iflag,iter,j,jm1,l,msum,ncfail,ncsuc,nslow1,nslow2
    integer iwa(1)
    logical jeval,sing
    double precision actred,delta,epsmch,fnorm,fnorm1,one,pnorm, &
         prered,p1,p5,p001,p0001,ratio,sum,temp,xnorm, &
         zero
!    double precision dpmpar,enorm
    data one,p1,p5,p001,p0001,zero &
         /1.0d0,1.0d-1,5.0d-1,1.0d-3,1.0d-4,0.0d0/
!
!     epsmch is the machine precision.
!
    epsmch = dpmpar(1)
!
    info = 0
    iflag = 0
    nfev = 0
!
!     check the input parameters for errors.
!
    if (n .le. 0 .or. xtol .lt. zero .or. maxfev .le. 0 &
         .or. ml .lt. 0 .or. mu .lt. 0 .or. factor .le. zero &
         .or. ldfjac .lt. n .or. lr .lt. (n*(n + 1))/2) go to 300
    if (mode .ne. 2) go to 20
    do  j = 1, n
       if (diag(j) .le. zero) go to 300
    enddo
20  continue
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
    iflag = 1
    call fcn(n,x,fvec,iflag)
    nfev = 1
    if (iflag .lt. 0) go to 300
    fnorm = enorm(n,fvec)
!
!     determine the number of calls to fcn needed to compute
!     the jacobian matrix.
!
    msum = min0(ml+mu+1,n)
!
!     initialize iteration counter and monitors.
!
    iter = 1
    ncsuc = 0
    ncfail = 0
    nslow1 = 0
    nslow2 = 0
!
!     beginning of the outer loop.
!
30  continue
    jeval = .true.
!
!        calculate the jacobian matrix.
!
    iflag = 2
    call fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1, &
         wa2)
    nfev = nfev + msum
    if (iflag .lt. 0) go to 300
!
!        compute the qr factorization of the jacobian.
!
    call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
    if (iter .ne. 1) go to 70
    if (mode .eq. 2) go to 50
    do j = 1, n
       diag(j) = wa2(j)
       if (wa2(j) .eq. zero) diag(j) = one
    enddo
50  continue
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
    do j = 1, n
       wa3(j) = diag(j)*x(j)
    enddo
    xnorm = enorm(n,wa3)
    delta = factor*xnorm
    if (delta .eq. zero) delta = factor
70  continue
!
!        form (q transpose)*fvec and store in qtf.
!
    do i = 1, n
       qtf(i) = fvec(i)
    enddo
    do j = 1, n
       if (fjac(j,j) .eq. zero) go to 110
       sum = zero
       do i = j, n
          sum = sum + fjac(i,j)*qtf(i)
       enddo
       temp = -sum/fjac(j,j)
       do i = j, n
          qtf(i) = qtf(i) + fjac(i,j)*temp
       enddo
110    continue
    enddo
!
!        copy the triangular factor of the qr factorization into r.
!
    sing = .false.
    do j = 1, n
       l = j
       jm1 = j - 1
       if (jm1 .lt. 1) go to 140
       do  i = 1, jm1
          r(l) = fjac(i,j)
          l = l + n - i
       enddo
140    continue
       r(l) = wa1(j)
       if (wa1(j) .eq. zero) sing = .true.
    enddo
!
!        accumulate the orthogonal factor in fjac.
!
    call qform(n,n,fjac,ldfjac,wa1)
!
!        rescale if necessary.
!
    if (mode .eq. 2) go to 170
    do j = 1, n
       diag(j) = dmax1(diag(j),wa2(j))
    enddo
170 continue
!
!        beginning of the inner loop.
!
180 continue
!
!           if requested, call fcn to enable printing of iterates.
!
    if (nprint .le. 0) go to 190
    iflag = 0
    if (mod(iter-1,nprint) .eq. 0) call fcn(n,x,fvec,iflag)
    if (iflag .lt. 0) go to 300
190 continue
!
!           determine the direction p.
!
    call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
!
!           store the direction p and x + p. calculate the norm of p.
!
    do j = 1, n
       wa1(j) = -wa1(j)
       wa2(j) = x(j) + wa1(j)
       wa3(j) = diag(j)*wa1(j)
    enddo
    pnorm = enorm(n,wa3)
!
!           on the first iteration, adjust the initial step bound.
!
    if (iter .eq. 1) delta = dmin1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
    iflag = 1
    call fcn(n,wa2,wa4,iflag)
    nfev = nfev + 1
    if (iflag .lt. 0) go to 300
    fnorm1 = enorm(n,wa4)
!
!           compute the scaled actual reduction.
!
    actred = -one
    if (fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction.
!
    l = 1
    do i = 1, n
       sum = zero
       do j = i, n
          sum = sum + r(l)*wa1(j)
          l = l + 1
       enddo
       wa3(i) = qtf(i) + sum
    enddo
    temp = enorm(n,wa3)
    prered = zero
    if (temp .lt. fnorm) prered = one - (temp/fnorm)**2
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
    ratio = zero
    if (prered .gt. zero) ratio = actred/prered
!
!           update the step bound.
!
    if (ratio .ge. p1) go to 230
    ncsuc = 0
    ncfail = ncfail + 1
    delta = p5*delta
    go to 240
230 continue
    ncfail = 0
    ncsuc = ncsuc + 1
    if (ratio .ge. p5 .or. ncsuc .gt. 1) &
         delta = dmax1(delta,pnorm/p5)
    if (dabs(ratio-one) .le. p1) delta = pnorm/p5
240 continue
!
!           test for successful iteration.
!
    if (ratio .lt. p0001) go to 260
!
!           successful iteration. update x, fvec, and their norms.
!
    do j = 1, n
       x(j) = wa2(j)
       wa2(j) = diag(j)*x(j)
       fvec(j) = wa4(j)
    enddo
    xnorm = enorm(n,wa2)
    fnorm = fnorm1
    iter = iter + 1
260 continue
!
!           determine the progress of the iteration.
!
    nslow1 = nslow1 + 1
    if (actred .ge. p001) nslow1 = 0
    if (jeval) nslow2 = nslow2 + 1
    if (actred .ge. p1) nslow2 = 0
!
!           test for convergence.
!
    if (delta .le. xtol*xnorm .or. fnorm .eq. zero) info = 1
    if (info .ne. 0) go to 300
!
!           tests for termination and stringent tolerances.
!
    if (nfev .ge. maxfev) info = 2
    if (p1*dmax1(p1*delta,pnorm) .le. epsmch*xnorm) info = 3
    if (nslow2 .eq. 5) info = 4
    if (nslow1 .eq. 10) info = 5
    if (info .ne. 0) go to 300
!
!           criterion for recalculating jacobian approximation
!           by forward differences.
!
    if (ncfail .eq. 2) go to 290
!
!           calculate the rank one modification to the jacobian
!           and update qtf if necessary.
!
    do j = 1, n
       sum = zero
       do i = 1, n
          sum = sum + fjac(i,j)*wa4(i)
       enddo
       wa2(j) = (sum - wa3(j))/pnorm
       wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
       if (ratio .ge. p0001) qtf(j) = sum
    enddo
!
!           compute the qr factorization of the updated jacobian.
!
    call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
    call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)
    call r1mpyq(1,n,qtf,1,wa2,wa3)
!
!           end of the inner loop.
!
    jeval = .false.
    go to 180
290 continue
!
!        end of the outer loop.
!
    go to 30
300 continue
!
!     termination, either normal or user imposed.
!
    if (iflag .lt. 0) info = iflag
    iflag = 0
    if (nprint .gt. 0) call fcn(n,x,fvec,iflag)
    return
!
!     last card of subroutine hybrd.
!
  end subroutine hybrd

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\


  subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
    implicit none
    integer n,lr
    double precision delta
    double precision r(lr),diag(n),qtb(n),x(n),wa1(n),wa2(n)
!     **********
!
!     subroutine dogleg
!
!     given an m by n matrix a, an n by n nonsingular diagonal
!     matrix d, an m-vector b, and a positive number delta, the
!     problem is to determine the convex combination x of the
!     gauss-newton and scaled gradient directions that minimizes
!     (a*x - b) in the least squares sense, subject to the
!     restriction that the euclidean norm of d*x be at most delta.
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization of a. that is, if a = q*r, where q has
!     orthogonal columns and r is an upper triangular matrix,
!     then dogleg expects the full upper triangle of r and
!     the first n components of (q transpose)*b.
!
!     the subroutine statement is
!
!       subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an input array of length lr which must contain the upper
!         triangular matrix r stored by rows.
!
!       lr is a positive integer input variable not less than
!         (n*(n+1))/2.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.
!
!       delta is a positive input variable which specifies an upper
!         bound on the euclidean norm of d*x.
!
!       x is an output array of length n which contains the desired
!         convex combination of the gauss-newton direction and the
!         scaled gradient direction.
!
!       wa1 and wa2 are work arrays of length n.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar,enorm
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
    integer i,j,jj,jp1,k,l
    double precision alpha,bnorm,epsmch,gnorm,one,qnorm,sgnorm,sum, &
         temp,zero
!    double precision dpmpar,enorm
    data one,zero /1.0d0,0.0d0/
!
!     epsmch is the machine precision.
!
    epsmch = dpmpar(1)
!
!     first, calculate the gauss-newton direction.
!
    jj = (n*(n + 1))/2 + 1
    do k = 1, n
       j = n - k + 1
       jp1 = j + 1
       jj = jj - k
       l = jj + 1
       sum = zero
       if (n .lt. jp1) go to 20
       do i = jp1, n
          sum = sum + r(l)*x(i)
          l = l + 1
       enddo
20     continue
       temp = r(jj)
       if (temp .ne. zero) go to 40
       l = j
       do i = 1, j
          temp = dmax1(temp,dabs(r(l)))
          l = l + n - i
       enddo
       temp = epsmch*temp
       if (temp .eq. zero) temp = epsmch
40     continue
       x(j) = (qtb(j) - sum)/temp
    enddo
!
!     test whether the gauss-newton direction is acceptable.
!
    do j = 1, n
       wa1(j) = zero
       wa2(j) = diag(j)*x(j)
    enddo
    qnorm = enorm(n,wa2)
    if (qnorm .le. delta) go to 140
!
!     the gauss-newton direction is not acceptable.
!     next, calculate the scaled gradient direction.
!
    l = 1
    do j = 1, n
       temp = qtb(j)
       do i = j, n
          wa1(i) = wa1(i) + r(l)*temp
          l = l + 1
       enddo
       wa1(j) = wa1(j)/diag(j)
    enddo
!
!     calculate the norm of the scaled gradient and test for
!     the special case in which the scaled gradient is zero.
!
    gnorm = enorm(n,wa1)
    sgnorm = zero
    alpha = delta/qnorm
    if (gnorm .eq. zero) go to 120
!
!     calculate the point along the scaled gradient
!     at which the quadratic is minimized.
!
    do j = 1, n
       wa1(j) = (wa1(j)/gnorm)/diag(j)
    enddo
    l = 1
    do j = 1, n
       sum = zero
       do i = j, n
          sum = sum + r(l)*wa1(i)
          l = l + 1
       enddo
       wa2(j) = sum
    enddo
    temp = enorm(n,wa2)
    sgnorm = (gnorm/temp)/temp
!
!     test whether the scaled gradient direction is acceptable.
!
    alpha = zero
    if (sgnorm .ge. delta) go to 120
!
!     the scaled gradient direction is not acceptable.
!     finally, calculate the point along the dogleg
!     at which the quadratic is minimized.
!
    bnorm = enorm(n,qtb)
    temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)
    temp = temp - (delta/qnorm)*(sgnorm/delta)**2 &
         + dsqrt((temp-(delta/qnorm))**2 &
         +(one-(delta/qnorm)**2)*(one-(sgnorm/delta)**2))
    alpha = ((delta/qnorm)*(one - (sgnorm/delta)**2))/temp
120 continue
!
!     form appropriate convex combination of the gauss-newton
!     direction and the scaled gradient direction.
!
    temp = (one - alpha)*dmin1(sgnorm,delta)
    do j = 1, n
       x(j) = temp*wa1(j) + alpha*x(j)
    enddo
140 continue
    return
!
!     last card of subroutine dogleg.
!
  end subroutine dogleg

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,&
       wa1,wa2)
    implicit none
    integer n,ldfjac,iflag,ml,mu
    double precision epsfcn
    double precision x(n),fvec(n),fjac(ldfjac,n),wa1(n),wa2(n)
! added external fcn
    external fcn
!     **********
!
!     subroutine fdjac1
!
!     this subroutine computes a forward-difference approximation
!     to the n by n jacobian matrix associated with a specified
!     problem of n functions in n variables. if the jacobian has
!     a banded form, then function evaluations are saved by only
!     approximating the nonzero terms.
!
!     the subroutine statement is
!
!       subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,
!                         wa1,wa2)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,iflag)
!         integer n,iflag
!         double precision x(n),fvec(n)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of fdjac1.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an input array of length n.
!
!       fvec is an input array of length n which must contain the
!         functions evaluated at x.
!
!       fjac is an output n by n array which contains the
!         approximation to the jacobian matrix evaluated at x.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       iflag is an integer variable which can be used to terminate
!         the execution of fdjac1. see description of fcn.
!
!       ml is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         ml to at least n - 1.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       mu is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         mu to at least n - 1.
!
!       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at
!         least n, then the jacobian is considered dense, and wa2 is
!         not referenced.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar
!
!       fortran-supplied ... dabs,dmax1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
    integer i,j,k,msum
    double precision eps,epsmch,h,temp,zero
!    double precision dpmpar
    data zero /0.0d0/
!
!     epsmch is the machine precision.
!
    epsmch = dpmpar(1)
!
    eps = dsqrt(dmax1(epsfcn,epsmch))
    msum = ml + mu + 1
    if (msum .lt. n) go to 40
!
!        computation of dense approximate jacobian.
!
    do j = 1, n
       temp = x(j)
       h = eps*dabs(temp)
       if (h .eq. zero) h = eps
       x(j) = temp + h
       call fcn(n,x,wa1,iflag)
       if (iflag .lt. 0) go to 30
       x(j) = temp
       do i = 1, n
          fjac(i,j) = (wa1(i) - fvec(i))/h
       enddo
    enddo
30  continue
    go to 110
40  continue
!
!        computation of banded approximate jacobian.
!
    do k = 1, msum
       do j = k, n, msum
          wa2(j) = x(j)
          h = eps*dabs(wa2(j))
          if (h .eq. zero) h = eps
          x(j) = wa2(j) + h
       enddo
       call fcn(n,x,wa1,iflag)
       if (iflag .lt. 0) go to 100
       do j = k, n, msum
          x(j) = wa2(j)
          h = eps*dabs(wa2(j))
          if (h .eq. zero) h = eps
          do i = 1, n
             fjac(i,j) = zero
             if (i .ge. j - mu .and. i .le. j + ml) &
                  fjac(i,j) = (wa1(i) - fvec(i))/h
          enddo
       enddo
    enddo
100 continue
110 continue
    return
!
!     last card of subroutine fdjac1.
!
  end subroutine fdjac1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  subroutine qform(m,n,q,ldq,wa)
    implicit none
    integer m,n,ldq
    double precision q(ldq,m),wa(m)
!     **********
!
!     subroutine qform
!
!     this subroutine proceeds from the computed qr factorization of
!     an m by n matrix a to accumulate the m by m orthogonal matrix
!     q from its factored form.
!
!     the subroutine statement is
!
!       subroutine qform(m,n,q,ldq,wa)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a and the order of q.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       q is an m by m array. on input the full lower trapezoid in
!         the first min(m,n) columns of q contains the factored form.
!         on output q has been accumulated into a square matrix.
!
!       ldq is a positive integer input variable not less than m
!         which specifies the leading dimension of the array q.
!
!       wa is a work array of length m.
!
!     subprograms called
!
!       fortran-supplied ... min0
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
    integer i,j,jm1,k,l,minmn,np1
    double precision one,sum,temp,zero
    data one,zero /1.0d0,0.0d0/
!
!     zero out upper triangle of q in the first min(m,n) columns.
!
    minmn = min0(m,n)
    if (minmn .lt. 2) go to 30
    do j = 2, minmn
       jm1 = j - 1
       do i = 1, jm1
          q(i,j) = zero
       enddo
    enddo
30  continue
!
!     initialize remaining columns to those of the identity matrix.
!
    np1 = n + 1
    if (m .lt. np1) go to 60
    do j = np1, m
       do i = 1, m
          q(i,j) = zero
       enddo
       q(j,j) = one
    enddo
60  continue
!
!     accumulate q from its factored form.
!
    do l = 1, minmn
       k = minmn - l + 1
       do i = k, m
          wa(i) = q(i,k)
          q(i,k) = zero
       enddo
       q(k,k) = one
       if (wa(k) .eq. zero) go to 110
       do j = k, m
          sum = zero
          do i = k, m
             sum = sum + q(i,j)*wa(i)
          enddo
          temp = sum/wa(k)
          do i = k, m
             q(i,j) = q(i,j) - temp*wa(i)
          enddo
       enddo
110    continue
    enddo
    return
!
!     last card of subroutine qform.
!
  end subroutine qform
      
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
    implicit none
    integer m,n,lda,lipvt
    integer ipvt(lipvt)
    logical pivot
    double precision a(lda,n),rdiag(n),acnorm(n),wa(n)
!     **********
!
!     subroutine qrfac
!
!     this subroutine uses householder transformations with column
!     pivoting (optional) to compute a qr factorization of the
!     m by n matrix a. that is, qrfac determines an orthogonal
!     matrix q, a permutation matrix p, and an upper trapezoidal
!     matrix r with diagonal elements of nonincreasing magnitude,
!     such that a*p = q*r. the householder transformation for
!     column k, k = 1,2,...,min(m,n), is of the form
!
!                           t
!           i - (1/u(k))*u*u
!
!     where u has zeros in the first k-1 positions. the form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding linpack subroutine.
!
!     the subroutine statement is
!
!       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       a is an m by n array. on input a contains the matrix for
!         which the qr factorization is to be computed. on output
!         the strict upper trapezoidal part of a contains the strict
!         upper trapezoidal part of r, and the lower trapezoidal
!         part of a contains a factored form of q (the non-trivial
!         elements of the u vectors described above).
!
!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.
!
!       pivot is a logical input variable. if pivot is set true,
!         then column pivoting is enforced. if pivot is set false,
!         then no column pivoting is done.
!
!       ipvt is an integer output array of length lipvt. ipvt
!         defines the permutation matrix p such that a*p = q*r.
!         column j of p is column ipvt(j) of the identity matrix.
!         if pivot is false, ipvt is not referenced.
!
!       lipvt is a positive integer input variable. if pivot is false,
!         then lipvt may be as small as 1. if pivot is true, then
!         lipvt must be at least n.
!
!       rdiag is an output array of length n which contains the
!         diagonal elements of r.
!
!       acnorm is an output array of length n which contains the
!         norms of the corresponding columns of the input matrix a.
!         if this information is not needed, then acnorm can coincide
!         with rdiag.
!
!       wa is a work array of length n. if pivot is false, then wa
!         can coincide with rdiag.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar,enorm
!
!       fortran-supplied ... dmax1,dsqrt,min0
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
    integer i,j,jp1,k,kmax,minmn
    double precision ajnorm,epsmch,one,p05,sum,temp,zero
!    double precision dpmpar,enorm
    data one,p05,zero /1.0d0,5.0d-2,0.0d0/
!
!     epsmch is the machine precision.
!
    epsmch = dpmpar(1)
!
!     compute the initial column norms and initialize several arrays.
!
    do j = 1, n
       acnorm(j) = enorm(m,a(1,j))
       rdiag(j) = acnorm(j)
       wa(j) = rdiag(j)
       if (pivot) ipvt(j) = j
    enddo
!10  continue
!
!     reduce a to r with householder transformations.
!
    minmn = min0(m,n)
!    do 110 j = 1, minmn
    do j = 1, minmn
       if (.not.pivot) go to 40
!
!        bring the column of largest norm into the pivot position.
!
       kmax = j
       do k = j, n
          if (rdiag(k) .gt. rdiag(kmax)) kmax = k
       enddo
!20     continue
       if (kmax .eq. j) go to 40
       do i = 1, m
          temp = a(i,j)
          a(i,j) = a(i,kmax)
          a(i,kmax) = temp
       enddo
!30     continue
       rdiag(kmax) = rdiag(j)
       wa(kmax) = wa(j)
       k = ipvt(j)
       ipvt(j) = ipvt(kmax)
       ipvt(kmax) = k
40     continue
!
!        compute the householder transformation to reduce the
!        j-th column of a to a multiple of the j-th unit vector.
!
       ajnorm = enorm(m-j+1,a(j,j))
       if (ajnorm .eq. zero) go to 100
       if (a(j,j) .lt. zero) ajnorm = -ajnorm
       do i = j, m
          a(i,j) = a(i,j)/ajnorm
       enddo
!50     continue
       a(j,j) = a(j,j) + one
!
!        apply the transformation to the remaining columns
!        and update the norms.
!
       jp1 = j + 1
       if (n .lt. jp1) go to 100
!       do 90 k = jp1, n
       do k = jp1, n
          sum = zero
          do i = j, m
             sum = sum + a(i,j)*a(i,k)
          enddo
!60        continue
          temp = sum/a(j,j)
          do i = j, m
             a(i,k) = a(i,k) - temp*a(i,j)
          enddo
!70        continue
          if (.not.pivot .or. rdiag(k) .eq. zero) go to 80
          temp = a(j,k)/rdiag(k)
          rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
          if (p05*(rdiag(k)/wa(k))**2 .gt. epsmch) go to 80
          rdiag(k) = enorm(m-j,a(jp1,k))
          wa(k) = rdiag(k)
80        continue
       enddo
!90     continue
100    continue
       rdiag(j) = -ajnorm
    enddo
!110 continue
    return
!
!     last card of subroutine qrfac.
!
  end subroutine qrfac

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  
  subroutine r1mpyq(m,n,a,lda,v,w)
    implicit none
    integer m,n,lda
    double precision a(lda,n),v(n),w(n)
!     **********
!
!     subroutine r1mpyq
!
!     given an m by n matrix a, this subroutine computes a*q where
!     q is the product of 2*(n - 1) transformations
!
!           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
!     and gv(i), gw(i) are givens rotations in the (i,n) plane which
!     eliminate elements in the i-th and n-th planes, respectively.
!     q itself is not given, rather the information to recover the
!     gv, gw rotations is supplied.
!
!     the subroutine statement is
!
!       subroutine r1mpyq(m,n,a,lda,v,w)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       a is an m by n array. on input a must contain the matrix
!         to be postmultiplied by the orthogonal matrix q
!         described above. on output a*q has replaced a.
!
!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.
!
!       v is an input array of length n. v(i) must contain the
!         information necessary to recover the givens rotation gv(i)
!         described above.
!
!       w is an input array of length n. w(i) must contain the
!         information necessary to recover the givens rotation gw(i)
!         described above.
!
!     subroutines called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
    integer i,j,nmj,nm1
    double precision cos,one,sin,temp
    data one /1.0d0/
!
!     apply the first set of givens rotations to a.
!
    nm1 = n - 1
    if (nm1 .lt. 1) go to 50
    do nmj = 1, nm1
       j = n - nmj
       if (dabs(v(j)) .gt. one) cos = one/v(j)
       if (dabs(v(j)) .gt. one) sin = dsqrt(one-cos**2)
       if (dabs(v(j)) .le. one) sin = v(j)
       if (dabs(v(j)) .le. one) cos = dsqrt(one-sin**2)
       do i = 1, m
          temp = cos*a(i,j) - sin*a(i,n)
          a(i,n) = sin*a(i,j) + cos*a(i,n)
          a(i,j) = temp
       enddo
    enddo
!
!     apply the second set of givens rotations to a.
!
    do j = 1, nm1
       if (dabs(w(j)) .gt. one) cos = one/w(j)
       if (dabs(w(j)) .gt. one) sin = dsqrt(one-cos**2)
       if (dabs(w(j)) .le. one) sin = w(j)
       if (dabs(w(j)) .le. one) cos = dsqrt(one-sin**2)
       do i = 1, m
          temp = cos*a(i,j) + sin*a(i,n)
          a(i,n) = -sin*a(i,j) + cos*a(i,n)
          a(i,j) = temp
       enddo
    enddo
50  continue
    return
!
!     last card of subroutine r1mpyq.
!
  end subroutine r1mpyq

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\


  subroutine r1updt(m,n,s,ls,u,v,w,sing)
    implicit none
    integer m,n,ls
    logical sing
    double precision s(ls),u(m),v(n),w(m)
!     **********
!
!     subroutine r1updt
!
!     given an m by n lower trapezoidal matrix s, an m-vector u,
!     and an n-vector v, the problem is to determine an
!     orthogonal matrix q such that
!
!                   t
!           (s + u*v )*q
!
!     is again lower trapezoidal.
!
!     this subroutine determines q as the product of 2*(n - 1)
!     transformations
!
!           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
!     where gv(i), gw(i) are givens rotations in the (i,n) plane
!     which eliminate elements in the i-th and n-th planes,
!     respectively. q itself is not accumulated, rather the
!     information to recover the gv, gw rotations is returned.
!
!     the subroutine statement is
!
!       subroutine r1updt(m,n,s,ls,u,v,w,sing)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of s.
!
!       n is a positive integer input variable set to the number
!         of columns of s. n must not exceed m.
!
!       s is an array of length ls. on input s must contain the lower
!         trapezoidal matrix s stored by columns. on output s contains
!         the lower trapezoidal matrix produced as described above.
!
!       ls is a positive integer input variable not less than
!         (n*(2*m-n+1))/2.
!
!       u is an input array of length m which must contain the
!         vector u.
!
!       v is an array of length n. on input v must contain the vector
!         v. on output v(i) contains the information necessary to
!         recover the givens rotation gv(i) described above.
!
!       w is an output array of length m. w(i) contains information
!         necessary to recover the givens rotation gw(i) described
!         above.
!
!       sing is a logical output variable. sing is set true if any
!         of the diagonal elements of the output s are zero. otherwise
!         sing is set false.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more,
!     john l. nazareth
!
!     **********
    integer i,j,jj,l,nmj,nm1
    double precision cos,cotan,giant,one,p5,p25,sin,tan,tau,temp, &
         zero
!    double precision dpmpar
    data one,p5,p25,zero /1.0d0,5.0d-1,2.5d-1,0.0d0/
!
!     giant is the largest magnitude.
!
    giant = dpmpar(3)
!
!     initialize the diagonal element pointer.
!
    jj = (n*(2*m - n + 1))/2 - (m - n)
!
!     move the nontrivial part of the last column of s into w.
!
    l = jj
    do i = n, m
       w(i) = s(l)
       l = l + 1
    enddo
!
!     rotate the vector v into a multiple of the n-th unit vector
!     in such a way that a spike is introduced into w.
!
    nm1 = n - 1
    if (nm1 .lt. 1) go to 70
    do nmj = 1, nm1
       j = n - nmj
       jj = jj - (m - j + 1)
       w(j) = zero
       if (v(j) .eq. zero) go to 50
!
!        determine a givens rotation which eliminates the
!        j-th element of v.
!
       if (dabs(v(n)) .ge. dabs(v(j))) go to 20
       cotan = v(n)/v(j)
       sin = p5/dsqrt(p25+p25*cotan**2)
       cos = sin*cotan
       tau = one
       if (dabs(cos)*giant .gt. one) tau = one/cos
       go to 30
20     continue
       tan = v(j)/v(n)
       cos = p5/dsqrt(p25+p25*tan**2)
       sin = cos*tan
       tau = sin
30     continue
!
!        apply the transformation to v and store the information
!        necessary to recover the givens rotation.
!
       v(n) = sin*v(j) + cos*v(n)
       v(j) = tau
!
!        apply the transformation to s and extend the spike in w.
!
       l = jj
       do i = j, m
          temp = cos*s(l) - sin*w(i)
          w(i) = sin*s(l) + cos*w(i)
          s(l) = temp
          l = l + 1
       enddo
50     continue
    enddo
70  continue
!
!     add the spike from the rank 1 update to w.
!
    do i = 1, m
       w(i) = w(i) + v(n)*u(i)
    enddo
!
!     eliminate the spike.
!
    sing = .false.
    if (nm1 .lt. 1) go to 140
    do j = 1, nm1
       if (w(j) .eq. zero) go to 120
!
!        determine a givens rotation which eliminates the
!        j-th element of the spike.
!
       if (dabs(s(jj)) .ge. dabs(w(j))) go to 90
       cotan = s(jj)/w(j)
       sin = p5/dsqrt(p25+p25*cotan**2)
       cos = sin*cotan
       tau = one
       if (dabs(cos)*giant .gt. one) tau = one/cos
       go to 100
90     continue
       tan = w(j)/s(jj)
       cos = p5/dsqrt(p25+p25*tan**2)
       sin = cos*tan
       tau = sin
100    continue
!
!        apply the transformation to s and reduce the spike in w.
!
       l = jj
       do i = j, m
          temp = cos*s(l) + sin*w(i)
          w(i) = -sin*s(l) + cos*w(i)
          s(l) = temp
          l = l + 1
       enddo
!
!        store the information necessary to recover the
!        givens rotation.
!
       w(j) = tau
120    continue
!
!        test for zero diagonal elements in the output s.
!
       if (s(jj) .eq. zero) sing = .true.
       jj = jj + (m - j + 1)
    enddo
140 continue
!
!     move w back into the last column of the output s.
!
    l = jj
    do i = n, m
       s(l) = w(i)
       l = l + 1
    enddo
    if (s(jj) .eq. zero) sing = .true.
    return
!
!     last card of subroutine r1updt.
!
  end subroutine r1updt

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  double precision function dpmpar(i)
    implicit none
    integer i
!     **********
!
!     Function dpmpar
!
!     This function provides double precision machine parameters
!     when the appropriate set of data statements is activated (by
!     removing the c from column 1) and all other data statements are
!     rendered inactive. Most of the parameter values were obtained
!     from the corresponding Bell Laboratories Port Library function.
!
!     The function statement is
!
!       double precision function dpmpar(i)
!
!     where
!
!       i is an integer input variable set to 1, 2, or 3 which
!         selects the desired machine parameter. If the machine has
!         t base b digits and its smallest and largest exponents are
!         emin and emax, respectively, then these parameters are
!
!         dpmpar(1) = b**(1 - t), the machine precision,
!
!         dpmpar(2) = b**(emin - 1), the smallest magnitude,
!
!         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
!
!     Argonne National Laboratory. MINPACK Project. November 1996.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'
!
!     **********
    integer mcheps(4)
    integer minmag(4)
    integer maxmag(4)
    double precision dmach(3)
    equivalence (dmach(1),mcheps(1))
    equivalence (dmach(2),minmag(1))
    equivalence (dmach(3),maxmag(1))
!
!     Machine constants for the IBM 360/370 series,
!     the Amdahl 470/V6, the ICL 2900, the Itel AS/6,
!     the Xerox Sigma 5/7/9 and the Sel systems 85/86.
!
!     data mcheps(1),mcheps(2) / z34100000, z00000000 /
!     data minmag(1),minmag(2) / z00100000, z00000000 /
!     data maxmag(1),maxmag(2) / z7fffffff, zffffffff /
!
!     Machine constants for the Honeywell 600/6000 series.
!
!     data mcheps(1),mcheps(2) / o606400000000, o000000000000 /
!     data minmag(1),minmag(2) / o402400000000, o000000000000 /
!     data maxmag(1),maxmag(2) / o376777777777, o777777777777 /
!
!     Machine constants for the CDC 6000/7000 series.
!
!     data mcheps(1) / 15614000000000000000b /
!     data mcheps(2) / 15010000000000000000b /
!
!     data minmag(1) / 00604000000000000000b /
!     data minmag(2) / 00000000000000000000b /
!
!     data maxmag(1) / 37767777777777777777b /
!     data maxmag(2) / 37167777777777777777b /
!
!     Machine constants for the PDP-10 (KA processor).
!
!     data mcheps(1),mcheps(2) / "114400000000, "000000000000 /
!     data minmag(1),minmag(2) / "033400000000, "000000000000 /
!     data maxmag(1),maxmag(2) / "377777777777, "344777777777 /
!
!     Machine constants for the PDP-10 (KI processor).
!
!     data mcheps(1),mcheps(2) / "104400000000, "000000000000 /
!     data minmag(1),minmag(2) / "000400000000, "000000000000 /
!     data maxmag(1),maxmag(2) / "377777777777, "377777777777 /
!
!     Machine constants for the PDP-11. 
!
!     data mcheps(1),mcheps(2) /   9472,      0 /
!     data mcheps(3),mcheps(4) /      0,      0 /
!
!     data minmag(1),minmag(2) /    128,      0 /
!     data minmag(3),minmag(4) /      0,      0 /
!
!     data maxmag(1),maxmag(2) /  32767,     -1 /
!     data maxmag(3),maxmag(4) /     -1,     -1 /
!
!     Machine constants for the Burroughs 6700/7700 systems.
!
!     data mcheps(1) / o1451000000000000 /
!     data mcheps(2) / o0000000000000000 /
!
!     data minmag(1) / o1771000000000000 /
!     data minmag(2) / o7770000000000000 /
!
!     data maxmag(1) / o0777777777777777 /
!     data maxmag(2) / o7777777777777777 /
!
!     Machine constants for the Burroughs 5700 system.
!
!     data mcheps(1) / o1451000000000000 /
!     data mcheps(2) / o0000000000000000 /
!
!     data minmag(1) / o1771000000000000 /
!     data minmag(2) / o0000000000000000 /
!
!     data maxmag(1) / o0777777777777777 /
!     data maxmag(2) / o0007777777777777 /
!
!     Machine constants for the Burroughs 1700 system.
!
!     data mcheps(1) / zcc6800000 /
!     data mcheps(2) / z000000000 /
!
!     data minmag(1) / zc00800000 /
!     data minmag(2) / z000000000 /
!
!     data maxmag(1) / zdffffffff /
!     data maxmag(2) / zfffffffff /
!
!     Machine constants for the Univac 1100 series.
!
!     data mcheps(1),mcheps(2) / o170640000000, o000000000000 /
!     data minmag(1),minmag(2) / o000040000000, o000000000000 /
!     data maxmag(1),maxmag(2) / o377777777777, o777777777777 /
!
!     Machine constants for the Data General Eclipse S/200.
!
!     Note - it may be appropriate to include the following card -
!     static dmach(3)
!
!     data minmag/20k,3*0/,maxmag/77777k,3*177777k/
!     data mcheps/32020k,3*0/
!
!     Machine constants for the Harris 220.
!
!     data mcheps(1),mcheps(2) / '20000000, '00000334 /
!     data minmag(1),minmag(2) / '20000000, '00000201 /
!     data maxmag(1),maxmag(2) / '37777777, '37777577 /
!
!     Machine constants for the Cray-1.
!
!     data mcheps(1) / 0376424000000000000000b /
!     data mcheps(2) / 0000000000000000000000b /
!
!     data minmag(1) / 0200034000000000000000b /
!     data minmag(2) / 0000000000000000000000b /
!
!     data maxmag(1) / 0577777777777777777777b /
!     data maxmag(2) / 0000007777777777777776b /
!
!     Machine constants for the Prime 400.
!
!     data mcheps(1),mcheps(2) / :10000000000, :00000000123 /
!     data minmag(1),minmag(2) / :10000000000, :00000100000 /
!     data maxmag(1),maxmag(2) / :17777777777, :37777677776 /
!
!     Machine constants for the VAX-11.
!
!     data mcheps(1),mcheps(2) /   9472,  0 /
!     data minmag(1),minmag(2) /    128,  0 /
!     data maxmag(1),maxmag(2) / -32769, -1 /
!
!     Machine constants for IEEE machines.
!
    data dmach(1) /2.22044604926d-16/
    data dmach(2) /2.22507385852d-308/
    data dmach(3) /1.79769313485d+308/
!
    dpmpar = dmach(i)
    return
!
!     Last card of function dpmpar.
!
  end function dpmpar

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  double precision function enorm(n,x)
    implicit none
    integer n
    double precision x(n)
!     **********
!
!     function enorm
!
!     given an n-vector x, this function calculates the
!     euclidean norm of x.
!
!     the euclidean norm is computed by accumulating the sum of
!     squares in three different sums. the sums of squares for the
!     small and large components are scaled so that no overflows
!     occur. non-destructive underflows are permitted. underflows
!     and overflows do not occur in the computation of the unscaled
!     sum of squares for the intermediate components.
!     the definitions of small, intermediate and large components
!     depend on two constants, rdwarf and rgiant. the main
!     restrictions on these constants are that rdwarf**2 not
!     underflow and rgiant**2 not overflow. the constants
!     given here are suitable for every known computer.
!
!     the function statement is
!
!       double precision function enorm(n,x)
!
!     where
!
!       n is a positive integer input variable.
!
!       x is an input array of length n.
!
!     subprograms called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
    integer i
    double precision agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs, &
         x1max,x3max,zero
    data one,zero,rdwarf,rgiant /1.0d0,0.0d0,3.834d-20,1.304d19/
    s1 = zero
    s2 = zero
    s3 = zero
    x1max = zero
    x3max = zero
    floatn = n
    agiant = rgiant/floatn
!    do 90 i = 1, n
    do i = 1, n
       xabs = dabs(x(i))
       if (xabs .gt. rdwarf .and. xabs .lt. agiant) go to 70
       if (xabs .le. rdwarf) go to 30
!
!              sum for large components.
!
       if (xabs .le. x1max) go to 10
       s1 = one + s1*(x1max/xabs)**2
       x1max = xabs
       go to 20
10     continue
       s1 = s1 + (xabs/x1max)**2
20     continue
       go to 60
30     continue
!
!              sum for small components.
!
       if (xabs .le. x3max) go to 40
       s3 = one + s3*(x3max/xabs)**2
       x3max = xabs
       go to 50
40     continue
       if (xabs .ne. zero) s3 = s3 + (xabs/x3max)**2
50     continue
60     continue
       go to 80
70     continue
!
!           sum for intermediate components.
!
       s2 = s2 + xabs**2
80     continue
    enddo
!90     continue
!
!     calculation of norm.
!
    if (s1 .eq. zero) go to 100
    enorm = x1max*dsqrt(s1+(s2/x1max)/x1max)
    go to 130
100 continue
    if (s2 .eq. zero) go to 110
    if (s2 .ge. x3max) &
         enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3))) 
    if (s2 .lt. x3max) &
         enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
    go to 120
110 continue
    enorm = x3max*dsqrt(s3)
120 continue
130 continue
    return
!
!     last card of function enorm.
!
  end function enorm

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

end MODULE MINPACK
