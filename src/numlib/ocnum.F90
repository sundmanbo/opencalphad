!
! MODULE LUKASNUM
MODULE OCNUM
!
! This is free software using LAPACK and BLAS
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
!
!-------------------
!
!  double precision, private, parameter :: DSMIN=1.0D-18,DSMAX=1.0D+18,&
!       DMAX=1.0D+36,DMIN=1.0D-36,epsx=1.0D-10,r=8.31451D0,&
!       unity=1.0D0,zero=0.0D0
! COMPILER WARNING ABOUT UNINITIALIZED ZERO, DSMIN, DSMAX, DMAX, DMIN, UNITY
!  double precision, private, parameter :: DSMIN=1.0D-18,DSMAX=1.0D+18,&
!       DMAX=1.0D+36,DMIN=1.0D-36,unity=1.0D0,zero=0.0D0
! modified to handle gas phases with fractions 1.0D-30 
!  double precision, private, parameter :: DSMIN=1.0D-33,DSMAX=1.0D+33,&
!       DMAX=1.0D+60,DMIN=1.0D-60,&
!       epsx=1.0D-10,r=8.31451D0,unity=1.0D0,zero=0.0D0
!
! The orignal routines LINGLD and MDINV written by H L Lukas
! has been replaced by using LAPACK and BLAS
!
! MDINV was split in two to handle symmetric and general matrices.
!
! oclablas is a small subset of lapack and blas needed for OC
! it is not optimized for any hardware.  If you have a full LAPACK+BLAS
! library for your hardware you should use that.
!
#ifdef NOLAPACK
use oclablas
! compile with -DLAPACK if LAPCK not extermal
#endif
!
! COMPILER WARNING ABOUT UNINITIALIZED ZERO, DSMIN, DSMAX, DMAX, DMIN, UNITY
  double precision, private, parameter :: DSMIN=1.0D-18,DSMAX=1.0D+18,&
       DMAX=1.0D+36,DMIN=1.0D-36,unity=1.0D0,zero=0.0D0
! declaration above must follow after USE
!
CONTAINS
    !
    !CCI
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    ! Development based on the work of Joao Pedro Carvalho Teuber 12/2020
    ! Linear system solved by splitting approch for conditions giving square
    ! mass matrix. Otherwise lingld is used
    SUBROUTINE lingldSplit(ND1,ND2,RMAT,X,N,M,NCONST,NPH)
        !-----------------------------------------------------------------------
        !     Solving a system of n linear equations with n unknowns
        !     | Masse  0      |  |X| = | Gibbs  |
        !     | Ca     MasseT |  | |   | Cig    |
        !-----------------------------------------------------------------------
        implicit none
        integer M,N,ND1,ND2, NPH, NCONST
        double precision RMAT(ND1,ND2),X(ND1)
        !-----------------------------------------------------------------------
        character trans*1
        integer i,j,k,nrhs,lda,ldb,info
        integer ipiv1(n), ipiv2(n), ipiv(n)
        double precision, allocatable :: a(:,:),Masse(:, :), Gibbs(:),MasseT(:, :), Cig(:)
        !
        allocate(a(n,n))
        allocate(Cig(NCONST))
        allocate(Masse(NPH, NCONST))
        allocate(Gibbs(NPH))
        allocate(MasseT(NCONST, NPH))
        !
        ipiv=0
        ipiv1=0
        ipiv2=0
        nrhs=1
        trans='N'
        lda=n
        ldb=n
        !
        do j=1,N
            do k=1,N
                a(j,k)=rmat(j,k)
            enddo
            x(j)=rmat(j,n+1)
        enddo
        do j=1,NPH
            do k=1, NCONST
                Masse(j,k) = a(j,k)
                Gibbs(j) = x(j)
            enddo
        enddo
        MasseT = transpose(Masse)
        ! Solve first part of the system
        call DGETRF(NPH, NCONST,Masse, NPH,IPIV1,INFO)
        if(info.ne.0) then
            write(*,*)'lingldSplit: Error return from dgetrf',info
            goto 900
        endif
        call DGETRS(TRANS,NPH,NRHS,Masse,NPH,IPIV1,Gibbs,NPH,INFO)
        ! Solve second part of the system
        ! Ca = a(j+NPH,1:NCONST)
        do j=1, NCONST
            Cig(j) = x(j+NPH) - DOT_PRODUCT(a(j+NPH,1:NCONST), Gibbs)
        enddo
        call DGETRF(NCONST, NPH, MasseT, NCONST,IPIV2,INFO)
        if(info.ne.0) then
            write(*,*)'lingldSplit: Error return from dgetrf',info
            goto 900
        endif
        call DGETRS(TRANS,NCONST,NRHS,MasseT,NCONST,IPIV2,Cig,NCONST,INFO)
        ! get solution
        do j=1, N
            if(j.LE.NCONST) then
                x(j) = Gibbs(j)
            else
                x(j) = Cig(j-NPH)
            endif
        enddo
900 continue
        deallocate(a,Masse, Gibbs,MasseT, Cig)
        return
    END SUBROUTINE lingldSplit
!CCI

!-----------------------------------------------------------------------
  SUBROUTINE LINGLDY (ND1,ND2,RMAT,X,N,M)
!-----------------------------------------------------------------------
!     System of n linear equations with n unknowns,
!     algorithm after Gauss with line exchange
!     ND1, ND2  =  Dimensioning of RMAT and X  (ND2 = ND1 + 1)
!     RMAT      =  matrix with right hand side as additional column, changed
!     X         =  result vector
!     N         =  number of equations and unknowns
!     M         =  Test for singularity (= n - rank)
!-----------------------------------------------------------------------
!---- COMMON VARIABLES
!    COMMON /ALLG/ DMAX,DMIN,DSMAX,DSMIN,EPSX,R,UNITY,ZERO
!    DOUBLE PRECISION DMAX,DMIN,DSMAX,DSMIN,EPSX,R,UNITY,ZERO
!-----------------------------------------------------------------------
!---- VARIABLES OF THE ARGUMENT LIST
    INTEGER M,N,ND1,ND2
    DOUBLE PRECISION RMAT(ND1,ND2),X(ND1)
!-----------------------------------------------------------------------
!---- LOCAL VARIABLES
    DOUBLE PRECISION A,B,C
    INTEGER I,I1,J,K,L,N1,NN1
!-----------------------------------------------------------------------
!    double precision, parameter :: DSMIN=1.0D-18,DSMAX=1.0D+18,&
!         DMAX=1.0D+36,DMIN=1.0D-36,epsx=1.0D-10,r=8.31451D0,&
!         unity=1.0D0,zero=0.0D0
!    DSMIN=1.0D-18
!    DSMAX=1.0D+18
!    DMAX=1.0D+36
!    DMIN=1.0D-36
!    epsx=1.0D-10
!    r=8.31451d0
!    unity=1.0D0
!    zero=0.0D0
!-----------------------------------------------------------------------
!    write(*,*)'enter lingld'
    N1=N+1
    NN1=N-1
    M=0
!-----------------------------------------------------------------------
    L490: DO I=1,NN1
       I1=I+1
       A=ZERO
       L=I
!-----------------------------------------------------------------------
!     search of pivot line
       L290: DO J=I,N
          B=ZERO
          L220: DO K=I,N
             IF(DABS(RMAT(J,K)).LT.DSMIN) GOTO 210
             IF(DABS(RMAT(J,K)).GT.DSMAX) GOTO 200
             B=B+RMAT(J,K)**2
             GOTO 220
200          B=B+DMAX
             GOTO 220
210          B=B+DMIN
220          CONTINUE
          enddo L220
          IF (B.GE.DMAX) GOTO 290
          IF (DABS(RMAT(J,I))*DSMAX.LT.B) GOTO 290
          IF (DABS(RMAT(J,I)).LT.DSMIN) GOTO 290
          C=RMAT(J,I)/B*RMAT(J,I)
          IF (C.LE.A) GOTO 290
          A=C
          L=J
290       CONTINUE
       enddo L290
!-----------------------------------------------------------------------
!     line exchange
       IF (L.EQ.I) GOTO 400
       L300: DO J=I,N1
          C=RMAT(I,J)
          RMAT(I,J)=RMAT(L,J)
          RMAT(L,J)=C
       enddo L300
!-----------------------------------------------------------------------
!     diagonalisation of the matrix
400    CONTINUE
       L470: DO J=I1,N
          IF (DABS(RMAT(J,I)).LT.DSMIN.AND.DABS(RMAT(I,I)).GE.UNITY) &
               GOTO 470
          IF (DABS(RMAT(J,I)).LE.UNITY.AND.DABS(RMAT(I,I)).GT.DSMAX) &
               GOTO 470
          IF (DABS(RMAT(J,I)).LT.DMIN) GOTO 470
          IF (DABS(RMAT(J,I)).GT.DSMAX.AND.DABS(RMAT(I,I)).LE.UNITY) &
               GOTO 460
          IF (DABS(RMAT(J,I)).GE.UNITY.AND.DABS(RMAT(I,I)).LT.DSMIN) &
               GOTO 460
          IF (DABS(RMAT(I,I)).LT.DMIN) GOTO 470
          C=RMAT(J,I)/RMAT(I,I)
          L440: DO K=I1,N1
             IF (DABS(RMAT(I,K)).LT.DSMIN.AND.DABS(C).LE.UNITY) GOTO 440
             IF (DABS(RMAT(I,K)).LE.UNITY.AND.DABS(C).LT.DSMIN) GOTO 440
             IF (DABS(RMAT(I,K)).LT.DMIN.OR.DABS(C).LT.DMIN) GOTO 440
             IF (DABS(RMAT(I,K)).GT.DSMAX.AND.DABS(C).GE.UNITY) GOTO 430
             IF (DABS(RMAT(I,K)).GE.UNITY.AND.DABS(C).GT.DSMAX) GOTO 430
             RMAT(J,K)=RMAT(J,K)-RMAT(I,K)*C
             GOTO 440
430          CALL WARNGB (3,A,B,J,K)
440          CONTINUE
          enddo L440
          GOTO 470
460       CALL WARNGB(3,A,B,J,I)
470       CONTINUE
       enddo L470
490    CONTINUE
    enddo L490
!-----------------------------------------------------------------------
    L700: DO L=1,N
       I=N1-L
       I1=I-1
       IF (DABS(RMAT(I,N1)).LT.DSMIN.AND.DABS(RMAT(I,I)).GE.UNITY) &
            GOTO 660
       IF (DABS(RMAT(I,N1)).LE.UNITY.AND.DABS(RMAT(I,I)).GT.DSMAX) &
            GOTO 660
       IF (DABS(RMAT(I,N1)).LT.DMIN) GOTO 660
       IF (DABS(RMAT(I,N1)).GT.DSMAX.AND.DABS(RMAT(I,I)).LE.UNITY) &
            GOTO 650
       IF (DABS(RMAT(I,N1)).GE.UNITY.AND.DABS(RMAT(I,I)).LT.DSMIN) &
            GOTO 650
       IF (DABS(RMAT(I,I)).LT.DMIN) GOTO 650
       C=RMAT(I,N1)/RMAT(I,I)
       X(I)=C
       IF (I.EQ.1) GOTO 700
       L600: DO J=1,I1
          IF (DABS(RMAT(J,I)).LT.DSMIN.AND.DABS(C).LE.UNITY) GOTO 600
          IF (DABS(RMAT(J,I)).LE.UNITY.AND.DABS(C).LT.DSMIN) GOTO 600
          IF (DABS(RMAT(J,I)).LT.DMIN.OR.DABS(C).LT.DMIN) GOTO 600
          IF (DABS(RMAT(J,I)).GT.DSMAX.AND.DABS(C).GE.UNITY) GOTO 550
          IF (DABS(RMAT(J,I)).GE.UNITY.AND.DABS(C).GT.DSMAX) GOTO 550
          RMAT(J,N1)=RMAT(J,N1)-RMAT(J,I)*C
          GOTO 600
550       CALL WARNGB (3,A,B,J,I)
600       CONTINUE
       enddo L600
       GOTO 700
!-----------------------------------------------------------------------
!     matrix singular, cutting of line and column to continue
650    M=M+1
660    X(I)=ZERO
!-----------------------------------------------------------------------
700    CONTINUE
    enddo L700
!    write(*,701)DSMIN,DSMAX,UNITY,DMAX,DMIN
!701 format('LINGLD: ',5(1PE12.4))
    RETURN
  END SUBROUTINE LINGLDY

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  SUBROUTINE WARNGB (NR,A1,A2,I1,I2)
!---- Printing of warnings, counting and stop printing after 5 times
!-----------------------------------------------------------------------
!---- VARIABLES OF THE ARGUMENT LIST
    INTEGER I1,I2,NR
    REAL*8 A1,A2
!-----------------------------------------------------------------------
!---- LOCAL VARIABLES
    INTEGER K(4)
    SAVE K
    DATA K/0,0,0,0/
!-----------------------------------------------------------------------
10  FORMAT (a,': Following message appears last time')
20  FORMAT (a,': Temperature',F8.2,' above maximum Temp.:',F8.2, &
    ' IPHEXC(*,1-2) =',2I3)
30  FORMAT (a,': Phase stability of component',I3,' of phase',I3/5X, &
    'is not defined for T =',F8.2,', range from',F8.2,' taken')
40  FORMAT (a,': Error in LINGLD, place(',I2,',',I2,')')
50  FORMAT (a,': d2G/dx2 suffers from rounding, phase',I3,' type',I3, &
    ' x =',E10.3,' test',E9.2)
90  FORMAT (a,': Subroutine WARNGB called with NR =',I3)
!-----------------------------------------------------------------------
!    write(*,*)'enter WARNGB'
    return
    IF (NR.LT.1.OR.NR.GT.4) GOTO 900
    K(NR)=K(NR)+1
    IF (K(NR).GT.5) RETURN
    IF (K(NR).EQ.5) WRITE (*,10)'WARNGB'
    GOTO (200,300,400,500),NR
200 WRITE (*,20)'WARNGB', A1,A2,I1,I2
    RETURN
300 WRITE (*,30)'WARNGB', I1,I2,A1,A2
    RETURN
400 WRITE (*,40)'WARNGB', I1,I2
    RETURN
500 WRITE (*,50)'WARNGB', I1,I2,A1,A2
    RETURN
900 WRITE (*,90)'WARNGB', NR
    STOP
  END SUBROUTINE WARNGB

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  SUBROUTINE PRECOND (ND1,ND2,RMAT,BADMAT)
! This is called from matsmin, code added by Clement/Joao from CEA
    implicit none
    INTEGER ND1,ND2, J
    DOUBLE PRECISION RMAT(ND1,ND2), INVERSIBLE
    LOGICAL BADMAT

    BADMAT=.FALSE.
    INVERSIBLE = 1.D20
    DO J=1,ND1
       INVERSIBLE = min(INVERSIBLE,DABS(RMAT(J,J)))
    ENDDO
    IF ((INVERSIBLE.GT.0.D0).AND.(ND1.EQ.ND2-1)) THEN
       DO J=1,ND1
          RMAT(J,ND2) = RMAT(J,ND2)/RMAT(J,J)
          RMAT(J,J) = 1.0D0
       ENDDO
    ELSE
! Added due to problems in parallel2 running all macros /2022.02.20 BOS
! probably because NEW Y does not reinitiate
!       write(*,*)'PRECOND: Matrix illconditioned',INVERSIBLE
       BADMAT=.TRUE.
! ignoring this message .... it does not seem to matter  2020.02.19/BoS
!        IF (ND1.NE.ND2-1) THEN
!            WRITE(*,*) 'PRECOND: No Square Matrix - no preconditiong applied'
!        ELSE
!           WRITE(*,77)ND1,ND2,INVERSIBLE
!77         format('PRECOND: Matrix not inversible - no preconditiong applied',&
!                2i4,1pe12.4)
!        ENDIF
    ENDIF
    RETURN
  END SUBROUTINE PRECOND

  SUBROUTINE LINGLD (ND1,ND2,RMAT,X,N,M)
!-----------------------------------------------------------------------
!     Solving a system of n linear equations with n unknowns
!     USING LAPACK+BLAS
!     ND1, ND2  =  Dimensioning of RMAT and X  (ND2 = ND1 + 1)
!     RMAT      =  matrix with right hand side as additional column, changed
!     X         =  result vector
!     N         =  number of equations and unknowns
!     M         =  Test for singularity (= n - rank)
!-----------------------------------------------------------------------
    implicit none
    INTEGER M,N,ND1,ND2
    DOUBLE PRECISION RMAT(ND1,ND2),X(ND1)
!-----------------------------------------------------------------------
    character trans*1
    integer j,k,nrhs,lda,ldb,info
    integer, allocatable :: ipiv(:)
    double precision, allocatable :: a(:,:)
!
    allocate(a(n,n))
    allocate(ipiv(n))
    ipiv=0
! there is just one right hand side
    nrhs=1
! right hand side is in rmat(n+1,j),j=1,n), move it to x
    do j=1,n
       do k=1,n
          a(j,k)=rmat(j,k)
       enddo
       x(j)=rmat(j,n+1)
    enddo
!    write(*,*)'Solving: ',nd1,nd2,n
!    do j=1,n
!       write(*,11)j,(rmat(j,k),k=1,n+1)
!    enddo
!    do j=1,n
!       write(*,11)j,x(j),(a(j,k),k=1,n)
!    enddo
11  format(i3,6(1pe12.4))
! trans='N' means no transpose
    trans='N'
    lda=n
    ldb=n
!
! we must first L*U factorize RMAT, the original values destroyed
!     CALL DGETRF(N,N,RMAT,LDA,IPIV,INFO)
     CALL DGETRF(N,N,A,LDA,IPIV,INFO)
     if(info.ne.0) then
!        write(*,*)'Error return dgetrf',info
        goto 900
     endif
! right hand side in X is overwritten by solution
     CALL DGETRS(TRANS,N,NRHS,A,LDA,IPIV,X,LDB,INFO)
!     if(info.ne.0) then
!        write(*,*)'Error return dgetrs',info
!     endif
900  continue
! info=0 meaks OK, returning m=0 means error
     m=info
! No warnings here, using gridminimizer may generate errors that can be ignored
!     if(info.gt.0) then
!        write(*,*)'Error solving equilibrium matrix with DGETRS'
!     else
!        write(*,*)'Solving equilibrium matrix with DGETRS',m
!     endif
1000 continue
!CCI
     deallocate (a, ipiv)
!CCI
     return
   END SUBROUTINE LINGLD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  SUBROUTINE MDINV (ND1,RMAT,RINV,N,IS)
! ND2 not used and removed to eliminate confusions
!  SUBROUTINE MDINV (ND1,ND2,RMAT,RINV,N,IS)
!  SUBROUTINE MSINV (ND1,RMAT,RINV,IS)
!-----------------------------------------------------------------------
!     Matrix inversion, symmetric matrix, DOUBLE PRECISION
!     using LAPACK (phase matrix)
!     ND1, ND2  =  Dimensioning of RMAT and RINV (ND2 = ND1 + 1)
!     RINV      =  invers of matrix RMAT (without last column of RMAT)
!     RMAT      =  matrix with additional column
!     N         =  number of lines and columns
!     IS        =  Test for singularity (0 = singular, 1 = not singular)
!-----------------------------------------------------------------------
!
    implicit none
    integer nd1,nd2,n,is
    double precision rmat(nd1,nd1),rinv(nd1,nd1)
!
    integer, dimension(:), allocatable :: ipiv
    double precision, dimension(:), allocatable :: work
    integer i,j,info,lda,m,lwork
    character uplo*1
!    if(nd1.ne.n) 
!    write(*,*)'in mdinv: ',nd1,n
! do not destroy RMAT
!    do i=1,nd1
!       write(*,17)nd1,(rmat(j,i),j=1,nd1)
!    enddo
!17  format(i3,6(1pe12.4)/(3x,6e12.4))
    RINV=RMAT
!
    lda=n
    allocate(ipiv(n))
! nonzero ipiv(i) will signalis the original row of row i
    ipiv=0
! upper triangular symmetric matrix
    uplo='U'
! if called with lwork=-1 the optimal dimension of work is returned
    m=-1
!    write(*,*)'Calling dsytrf',lda,m,n
    allocate(work(800))
    CALL DSYTRF(UPLO,N,RMAT,LDA,IPIV,WORK,m,INFO)
    if(info.ne.0) then
       write(*,*)'MDINV: Error from DSYTRF: ',info
       IS=0
       goto 1000
    endif
!
    lwork=int(work(1))
!    write(*,*)'lwork: ',nd1,n,lwork
    if(lwork.gt.700) then
       deallocate(work)
       allocate(work(lwork))
    endif
! factorize a symmetric unpacked indefinite matrix
    CALL DSYTRF(UPLO,N,RINV,LDA,IPIV,WORK,LWORK,INFO)
    if(info.ne.0) then
!       write(*,*)'Error return from DSYTRF:',info
       is=0; goto 1000
    endif
! invert using the factorization
    CALL DSYTRI(UPLO,N,RINV,LDA,IPIV,WORK,INFO)
!    write(*,*)'Info: ',info,n,lda,lwork
    if(info.ne.0) then
!       write(*,*)'Error return from DSYTRI: ',info
       is=0; goto 1000
    endif
! copy solution to RINV triangle to lower
  do i=2,n
     do j=1,i-1
        RINV(i,j)=RINV(j,i)
     enddo
  enddo
! all OK
!    write(*,*)'Matrix inverted using DSYTRI'
    is=1
!
1000 continue
    RETURN
  END SUBROUTINE MDINV

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  SUBROUTINE MDINVOLD(ND1,RMAT,RINV,N,IS)
! ND2 not used and removed to eliminate confusion
!  SUBROUTINE MDINVOLD(ND1,ND2,RMAT,RINV,N,IS)
!  SUBROUTINE MGINV(ND1,RMAT,RINV,IS)
!-----------------------------------------------------------------------
!     Matrix inversion, general matrix, DOUBLE PRECISION
!     using LAPACK for general matrix (component matrix)
!     ND1, ND2  =  Dimensioning of RMAT and RINV (ND2 = ND1 + 1)
!     RINV      =  invers of matrix RMAT (without last column of RMAT)
!     RMAT      =  matrix with additional column
!     N         =  number of lines and columns
!     IS        =  Test for singularity (0 = singular, 1 = not singular)
!-----------------------------------------------------------------------
!
    implicit none
    integer nd1,nd2,n,is
    double precision rmat(nd1,nd1),rinv(nd1,nd1)
    integer, dimension(:), allocatable :: ipiv
    double precision, dimension(:), allocatable :: work
    integer i,j,info,lda,m,lwork
!    if(nd1.ne.n) 
!    write(*,*)'in mdinv: ',nd1,n
! do not destroy RMAT
!    do i=1,nd1
!       write(*,17)nd1,(rmat(j,i),j=1,nd1)
!    enddo
!17  format(i3,6(1pe12.4)/(3x,6e12.4))
! copy input matrix to solution not to destroy RMAT
    RINV=RMAT
!
    lda=n
    allocate(ipiv(n))
! nonzero ipiv(i) will signal the original row of row i
    ipiv=0
! if called with lwork=-1 the optimal dimension of work is returned
    m=-1
!    write(*,*)'Calling dsytrf',lda,m,n
    allocate(work(800))
! replaced DSY with DGE for general matrix inversion .... ????
    CALL DGETRI(N,RINV,LDA,IPIV,WORK,m,INFO)
    if(info.ne.0) then
!       write(*,*)'Error from DGETRI at 1: ',info
       IS=0
       goto 1000
    endif
!
    lwork=int(work(1))
!    write(*,*)'lwork: ',nd1,n,lwork
    if(lwork.gt.700) then
       deallocate(work)
       allocate(work(lwork))
    endif
! factorize an general matrix
    CALL DGETRF(N,N,RINV,LDA,IPIV,INFO)
    if(info.ne.0) then
!       write(*,*)'Error return from DGETRF:',info
       is=0; goto 1000
    endif
! invert a general matrix using the factorization
    CALL DGETRI(N,RINV,LDA,IPIV,WORK,LWORK,INFO)
!    write(*,*)'Info: ',info,n,lda,lwork
    if(info.ne.0) then
!       write(*,*)'Error return from DGETRI: ',info
       is=0; goto 1000
    endif
!    write(*,*)'Matrix inverted using DGETRI'
! All OK, the solution is in RINV
    is=1
!
1000 continue
    return
  end SUBROUTINE MDINVOLD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

END MODULE OCNUM
