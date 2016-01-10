!
MODULE LUKASNUM
!
! This is free software developed by H L Lukas, Germany
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
!--------------------------------------------------------------------------
!
!  double precision, private, parameter :: DSMIN=1.0D-18,DSMAX=1.0D+18,&
!       DMAX=1.0D+36,DMIN=1.0D-36,epsx=1.0D-10,r=8.31451D0,&
!       unity=1.0D0,zero=0.0D0
! modified to handle gas phases with fractions 1.0D-30 
  double precision, private, parameter :: DSMIN=1.0D-33,DSMAX=1.0D+33,&
       DMAX=1.0D+60,DMIN=1.0D-60,&
       epsx=1.0D-10,r=8.31451D0,unity=1.0D0,zero=0.0D0
!  double precision, private :: DSMAX=1.0D+18
!  double precision, private :: DMAX=1.0D+36
!  double precision, private :: DMIN=1.0D-36
!  double precision, private :: epsx=1.0D-10
!  double precision, private :: r=8.31451d0
!  double precision, private :: unity=1.0D0
!  double precision, private :: zero=0.0D0
!
CONTAINS
!
  SUBROUTINE LINGLD (ND1,ND2,RMAT,X,N,M)
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
460       CALL WARNGB (3,A,B,J,I)
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
  END SUBROUTINE LINGLD

!///////////////////////////////////////////////////////////////////////

  SUBROUTINE MDINV (ND1,ND2,RMAT,RINV,N,IS)
!-----------------------------------------------------------------------
!     Matrix inversion, DOUBLE PRECISION
!     algorithm after Gauss for linear equations with line exchange
!     ND1, ND2  =  Dimensioning of RMAT and RINV (ND2 = ND1 + 1)
!     RINV      =  invers of matrix RMAT (without last column of RMAT)
!     RMAT      =  matrix with additional column
!     N         =  number of lines and columns
!     IS        =  Test for singularity (0 = singular, 1 = not singular)
!-----------------------------------------------------------------------
!---- COMMON VARIABLES
!    COMMON /ALLG/ DMAX,DMIN,DSMAX,DSMIN,EPSX,R,UNITY,ZERO
!    DOUBLE PRECISION DMAX,DMIN,DSMAX,DSMIN,EPSX,R,UNITY,ZERO
!-----------------------------------------------------------------------
!---- VARIABLES OF THE ARGUMENT LIST
    INTEGER N,IS,ND1,ND2
    DOUBLE PRECISION RMAT(ND1,ND2),RINV(ND1,ND1)
!-----------------------------------------------------------------------
!---- LOCAL VARIABLES
    DOUBLE PRECISION A,B,C
    INTEGER I,I1,J,K,L,N1
!    double precision, parameter :: DSMIN=1.0D-18,DSMAX=1.0D+18,&
!         DMAX=1.0D+36,DMIN=1.0D-36,epsx=1.0D-10,r=8.31451D0,&
!         unity=1.0D0,zero=0.0D0
!-----------------------------------------------------------------------
!    write(*,*)'enter mdinv'
    IS=1
    N1=N-1
    L110: DO I=1,N
       L100: DO J=1,N
          RINV(I,J)=ZERO
       enddo L100
       RINV(I,I)=UNITY
    enddo L110
    IF (N.LE.1) GOTO 500
!-----------------------------------------------------------------------
    L490: DO I=1,N1
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
          IF (B.GE.DMAX) GOTO 700
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
       IF (L.EQ.I) GOTO 330
       L310: DO J=I,N
          C=RMAT(I,J)
          RMAT(I,J)=RMAT(L,J)
          RMAT(L,J)=C
       enddo L310
       L320: DO J=1,N
          C=RINV(I,J)
          RINV(I,J)=RINV(L,J)
          RINV(L,J)=C
       enddo L320
!-----------------------------------------------------------------------
!     diagonalisation of the matrix
330    CONTINUE
       L480: DO J=I1,N
          IF (DABS(RMAT(J,I)).LT.DSMIN.AND.DABS(RMAT(I,I)).GE.UNITY) &
               GOTO 480
          IF (DABS(RMAT(J,I)).LE.UNITY.AND.DABS(RMAT(I,I)).GT.DSMAX) &
               GOTO 480
          IF (DABS(RMAT(J,I)).LT.DMIN) GOTO 480
          IF (DABS(RMAT(J,I)).GT.DSMAX.AND.DABS(RMAT(I,I)).LE.UNITY) &
               GOTO 700
          IF (DABS(RMAT(J,I)).GE.UNITY.AND.DABS(RMAT(I,I)).LT.DSMIN) &
               GOTO 700
          IF (DABS(RMAT(I,I)).LT.DMIN) GOTO 700
          C=RMAT(J,I)/RMAT(I,I)
          L440: DO K=I1,N
             IF (DABS(RMAT(I,K)).LT.DSMIN.AND.DABS(C).LE.UNITY) GOTO 440
             IF (DABS(RMAT(I,K)).LE.UNITY.AND.DABS(C).LT.DSMIN) GOTO 440
             IF (DABS(RMAT(I,K)).LT.DMIN.OR.DABS(C).LT.DMIN) GOTO 440
             IF (DABS(RMAT(I,K)).GT.DSMAX.AND.DABS(C).GE.UNITY) GOTO 700
             IF (DABS(RMAT(I,K)).GE.UNITY.AND.DABS(C).GT.DSMAX) GOTO 700
             RMAT(J,K)=RMAT(J,K)-RMAT(I,K)*C
440          CONTINUE
          enddo L440
          L460: DO K=1,N
             IF (DABS(RINV(I,K)).LT.DSMIN.AND.DABS(C).LE.UNITY) GOTO 460
             IF (DABS(RINV(I,K)).LE.UNITY.AND.DABS(C).LT.DSMIN) GOTO 460
             IF (DABS(RINV(I,K)).LT.DMIN.OR.DABS(C).LT.DMIN) GOTO 460
             IF (DABS(RINV(I,K)).GT.DSMAX.AND.DABS(C).GE.UNITY) GOTO 700
             IF (DABS(RINV(I,K)).GE.UNITY.AND.DABS(C).GT.DSMAX) GOTO 700
             RINV(J,K)=RINV(J,K)-RINV(I,K)*C
460          CONTINUE
          enddo L460
480       CONTINUE
       enddo L480
490    CONTINUE
    enddo L490
!-----------------------------------------------------------------------
500 CONTINUE
    L690A: DO K=1,N
       L690B: DO L=1,N
          I=N+1-L
          I1=I-1
          IF (DABS(RINV(I,K)).LT.DSMIN.AND.DABS(RMAT(I,I)).GE.UNITY) &
               GOTO 660
          IF (DABS(RINV(I,K)).LE.UNITY.AND.DABS(RMAT(I,I)).GT.DSMAX) &
               GOTO 660
          IF (DABS(RINV(I,K)).LT.DMIN) GOTO 660
          IF (DABS(RINV(I,K)).GT.DSMAX.AND.DABS(RMAT(I,I)).LE.UNITY) &
               GOTO 700
          IF (DABS(RINV(I,K)).GE.UNITY.AND.DABS(RMAT(I,I)).LT.DSMIN) &
               GOTO 700
          IF (DABS(RMAT(I,I)).LT.DMIN) GOTO 700
          RINV(I,K)=RINV(I,K)/RMAT(I,I)
          C=RINV(I,K)
          IF (I.EQ.1) GOTO 690
          L600: DO J=1,I1
             IF (DABS(RMAT(J,I)).LT.DSMIN.AND.DABS(C).LE.UNITY) GOTO 600
             IF (DABS(RMAT(J,I)).LE.UNITY.AND.DABS(C).LT.DSMIN) GOTO 600
             IF (DABS(RMAT(J,I)).LT.DMIN.OR.DABS(C).LT.DMIN) GOTO 600
             IF (DABS(RMAT(J,I)).GT.DSMAX.AND.DABS(C).GE.UNITY) GOTO 700
             IF (DABS(RMAT(J,I)).GE.UNITY.AND.DABS(C).GT.DSMAX) GOTO 700
             RINV(J,K)=RINV(J,K)-RMAT(J,I)*C
600          CONTINUE
          enddo L600
          GOTO 690
660       RINV(I,K)=ZERO
690       CONTINUE
       enddo L690B
    enddo L690A
    RETURN
!-----------------------------------------------------------------------
!     overflow ( = matrix singular)
700 IS=0
!    write(*,701)DSMIN,DSMAX,UNITY,DMAX,DMIN
!701 format('MDINV: ',5(1PE12.4))
    RETURN
  END SUBROUTINE MDINV
!///////////////////////////////////////////////////////////////////////
  SUBROUTINE DBLDET (HMAT,NDIM,N,DET)
!-----------------------------------------------------------------------
!   CALCULATING DETERMINANT OF MATRIX HMAT, DOUBLE PRECISION
!   BY GAUSS ALGORITHM FOR LINEAR EQUATIONS WITH LINE EXCHANGE
!   N   = NUMBER OF LINES AND COLUMNS OF MATRIX
!   DET = VALUE OF DETERMINANT
!   H.L. LUKAS, MAX-PLANCK-INSTITUT F. METALLFORSCHUNG,
!   HEISENBERGSTR. 5, D-7000 STUTTGART 80, F.R. GERMANY
!   LAST UPDATE: 10 SEPTEMBER 1986
!-----------------------------------------------------------------------
    DOUBLE PRECISION HMAT(NDIM,NDIM),A,B,C,DET
!    write(*,*)'enter dbldet'
    DET=1.0D+00
!-----------------------------------------------------------------------
    L490: DO I=1,N
       I1=I+1
       A=0.0D+00
       L=I
       IF (I.EQ.N) GOTO 300
!-----------------------------------------------------------------------
!---- SEARCH OF PIVOT LINE
       L290: DO J=I,N
          B=0.0D+00
          L220: DO K=I,N
             IF(DABS(HMAT(J,K)).LT.1.0D-36) GOTO 210
             IF(DABS(HMAT(J,K)).GT.1.0D+36) GOTO 200
             B=B+HMAT(J,K)**2
             GOTO 220
200          B=B+1.0D+72
             GOTO 220
210          B=B+1.0D-72
220          CONTINUE
          enddo L220
          IF (B.GE.1.0D+72) GOTO 520
          IF (DABS(HMAT(J,I))*1.0D+36.LT.B) GOTO 290
          C=HMAT(J,I)/B*HMAT(J,I)
          IF (C.LE.A) GOTO 290
          A=C
          L=J
290       CONTINUE
       enddo L290
!-----------------------------------------------------------------------
!---- LINE EXCHANGE
300    IF (L.EQ.I) GOTO 330
       L310: DO J=I,N
          C=HMAT(I,J)
          HMAT(I,J)=HMAT(L,J)
          HMAT(L,J)=C
       enddo L310
       DET=-DET
330    IF (DABS(DET).LT.1.0D-36.AND.DABS(HMAT(I,I)).LT.1.0D+00) GOTO 510
       IF (DABS(DET).LT.1.0D+00.AND.DABS(HMAT(I,I)).LT.1.0D-36) GOTO 510
       IF (DABS(DET).GT.1.0D+36.AND.DABS(HMAT(I,I)).GT.1.0D+00) GOTO 520
       IF (DABS(DET).GT.1.0D+00.AND.DABS(HMAT(I,I)).GT.1.0D+36) GOTO 520
       DET=DET*HMAT(I,I)
!-----------------------------------------------------------------------
!---- DIAGONALISATION OF THE MATRIX
400    IF (I.GE.N) GOTO 490
       L480: DO J=I1,N
          IF (DABS(HMAT(J,I)).LT.1.0D-36.AND.DABS(HMAT(I,I)).GT.1.0D+00) &
               GOTO 480
          IF (DABS(HMAT(J,I)).LT.1.0D+00.AND.DABS(HMAT(I,I)).GT.1.0D+36) &
               GOTO 480
          IF (DABS(HMAT(J,I)).GT.1.0D+36.AND.DABS(HMAT(I,I)).LT.1.0D+00) &
               GOTO 520
          IF (DABS(HMAT(J,I)).GT.1.0D+00.AND.DABS(HMAT(I,I)).LT.1.0D-36) &
               GOTO 520
          IF (DABS(HMAT(I,I)).LT.1.0D-72) GOTO 520
          C=HMAT(J,I)/HMAT(I,I)
          L440: DO K=I1,N
             IF (DABS(HMAT(I,K)).LT.1.0D-36.AND.DABS(C).LT.1.0D+00) GOTO 440
             IF (DABS(HMAT(I,K)).LT.1.0D+00.AND.DABS(C).LT.1.0D-36) GOTO 440
             IF (DABS(HMAT(I,K)).GT.1.0D+36.AND.DABS(C).GT.1.0D+00) GOTO 520
             IF (DABS(HMAT(I,K)).GT.1.0D+00.AND.DABS(C).GT.1.0D+36) GOTO 520
             HMAT(J,K)=HMAT(J,K)-HMAT(I,K)*C
440          CONTINUE
          enddo L440
480       CONTINUE
       enddo L480
490    CONTINUE
    enddo L490
!-----------------------------------------------------------------------
!---- CALCULATION SUCCESSFUL
500 RETURN
!-----------------------------------------------------------------------
!---- MATRIX SINGULAR
510 DET=0.0D+00
    RETURN
!-----------------------------------------------------------------------
!---- CALCULATION NOT POSSIBLE BECAUSE OF OVERFLOW
520 DET=1.0D+73
    RETURN
  END SUBROUTINE DBLDET
!/////////////////////////////////////////////////////////////////////
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
10  FORMAT (' Following message appears last time')
20  FORMAT (' Temperature',F8.2,' above maximum Temp.:',F8.2, &
    ' IPHEXC(*,1-2) =',2I3)
30  FORMAT (' Phase stability of component',I3,' of phase',I3/5X, &
    'is not defined for T =',F8.2,', range from',F8.2,' taken')
40  FORMAT (' Error in LINGLD, place(',I2,',',I2,')')
50  FORMAT (' d2G/dx2 suffers from rounding, phase',I3,' type',I3, &
    ' x =',E10.3,' test',E9.2)
90  FORMAT (' Subroutine WARNGB called with NR =',I3)
!-----------------------------------------------------------------------
!    write(*,*)'enter warngb'
    return
    IF (NR.LT.1.OR.NR.GT.4) GOTO 900
    K(NR)=K(NR)+1
    IF (K(NR).GT.5) RETURN
    IF (K(NR).EQ.5) WRITE (*,10)
    GOTO (200,300,400,500),NR
200 WRITE (*,20) A1,A2,I1,I2
    RETURN
300 WRITE (*,30) I1,I2,A1,A2
    RETURN
400 WRITE (*,40) I1,I2
    RETURN
500 WRITE (*,50) I1,I2,A1,A2
    RETURN
900 WRITE (*,90) NR
    STOP
  END SUBROUTINE WARNGB
!///////////////////////////////////////////////////////////////////////
!
END MODULE LUKASNUM
