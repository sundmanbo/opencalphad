!
MODULE OCLABLAS
!
! This is an extract of a fews routines from LAPACK and BLAS version 3.6.0
! used to invert a symmetric matrix and to solve a system of linear equations
! and some more things in DOUBLE PRECISION used in OpenCalphad
!
! LAPACk and BLAS are free software libraries
! Both converted from F77 to F90 in a minimal way (comments and 
! continuation lines modified).
!
CONTAINS
!
! -------------------------------------------------------------------------
!
! Original LAPACK/BLAS routines below
!
!
!> \brief \b DGETRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DGETRI + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetri.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetri.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetri.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGETRI computes the inverse of a matrix using the LU factorization
!> computed by DGETRF.
!>
!> This method inverts U and then computes inv(A) by solving the system
!> inv(A)*L = inv(U) for inv(A).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the factors L and U from the factorization
!>          A = P*L*U as computed by DGETRF.
!>          On exit, if INFO = 0, the inverse of the original matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,N).
!>          For optimal performance LWORK >= N*NB, where NB is
!>          the optimal blocksize returned by ILAENV.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!>                singular and its inverse could not be computed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup doubleGEcomputational
!
!  =====================================================================
  SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
    INTEGER            IPIV( * )
    DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
    DOUBLE PRECISION   ZERO, ONE
    PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
    LOGICAL            LQUERY
    INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB,&
         NBMIN, NN
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ILAENV
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEMM, DGEMV, DSWAP, DTRSM, DTRTRI, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
    INFO = 0
    NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
    LWKOPT = N*NB
    WORK( 1 ) = LWKOPT
    LQUERY = ( LWORK.EQ.-1 )
    IF( N.LT.0 ) THEN
       INFO = -1
    ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
       INFO = -3
    ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
       INFO = -6
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DGETRI', -INFO )
       RETURN
    ELSE IF( LQUERY ) THEN
       RETURN
    END IF
!
!     Quick return if possible
!
    IF( N.EQ.0 ) RETURN
!
!     Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
!     and the inverse is not computed.
!
    CALL DTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
    IF( INFO.GT.0 ) RETURN
!
    NBMIN = 2
    LDWORK = N
    IF( NB.GT.1 .AND. NB.LT.N ) THEN
       IWS = MAX( LDWORK*NB, 1 )
       IF( LWORK.LT.IWS ) THEN
          NB = LWORK / LDWORK
          NBMIN = MAX( 2, ILAENV( 2, 'DGETRI', ' ', N, -1, -1, -1 ) )
       END IF
    ELSE
       IWS = N
    END IF
!
!     Solve the equation inv(A)*L = inv(U) for inv(A).
!
    IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
!
!        Use unblocked code.
!
       DO 20 J = N, 1, -1
!
!           Copy current column of L to WORK and replace with zeros.
!
          DO 10 I = J + 1, N
             WORK( I ) = A( I, J )
             A( I, J ) = ZERO
10        CONTINUE
!
!           Compute current column of inv(A).
!
          IF( J.LT.N ) &
               CALL DGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),&
               LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
20     CONTINUE
    ELSE
!
!        Use blocked code.
!
       NN = ( ( N-1 ) / NB )*NB + 1
       DO 50 J = NN, 1, -NB
          JB = MIN( NB, N-J+1 )
!
!           Copy current block column of L to WORK and replace with
!           zeros.
!
          DO 40 JJ = J, J + JB - 1
             DO 30 I = JJ + 1, N
                WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                A( I, JJ ) = ZERO
30           CONTINUE
40        CONTINUE
!
!           Compute current block column of inv(A).
!
          IF( J+JB.LE.N ) &
               CALL DGEMM( 'No transpose', 'No transpose', N, JB,&
               N-J-JB+1, -ONE, A( 1, J+JB ), LDA,&
               WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
          CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,&
               ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
50     CONTINUE
    END IF
!
!     Apply column interchanges.
!
    DO 60 J = N - 1, 1, -1
       JP = IPIV( J )
       IF( JP.NE.J ) CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
60  CONTINUE
!
    WORK( 1 ) = IWS
    RETURN
!
!     End of DGETRI
!
 END SUBROUTINE DGETRI
!
!=
!
!> \brief \b DTRTRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DTRTRI + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrtri.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrtri.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrtri.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTRTRI computes the inverse of a real upper or lower triangular
!> matrix A.
!>
!> This is the Level 3 BLAS version of the algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  A is upper triangular;
!>          = 'L':  A is lower triangular.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          = 'N':  A is non-unit triangular;
!>          = 'U':  A is unit triangular.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the triangular matrix A.  If UPLO = 'U', the
!>          leading N-by-N upper triangular part of the array A contains
!>          the upper triangular matrix, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of the array A contains
!>          the lower triangular matrix, and the strictly upper
!>          triangular part of A is not referenced.  If DIAG = 'U', the
!>          diagonal elements of A are also not referenced and are
!>          assumed to be 1.
!>          On exit, the (triangular) inverse of the original matrix, in
!>          the same storage format.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
!>               matrix is singular and its inverse can not be computed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
 SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
   CHARACTER          DIAG, UPLO
   INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE, ZERO
   PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            NOUNIT, UPPER
   INTEGER            J, JB, NB, NN
!     ..
!     .. External Functions ..
!   LOGICAL            LSAME
!   INTEGER            ILAENV
!   EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
!   EXTERNAL           DTRMM, DTRSM, DTRTI2, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   NOUNIT = LSAME( DIAG, 'N' )
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
      INFO = -2
   ELSE IF( N.LT.0 ) THEN
      INFO = -3
   ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
      INFO = -5
   END IF
   IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'DTRTRI', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N.EQ.0 ) RETURN
!
!     Check for singularity if non-unit.
!
   IF( NOUNIT ) THEN
      DO 10 INFO = 1, N
         IF( A( INFO, INFO ).EQ.ZERO ) RETURN
10    CONTINUE
      INFO = 0
   END IF
!
!     Determine the block size for this environment.
!
   NB = ILAENV( 1, 'DTRTRI', UPLO // DIAG, N, -1, -1, -1 )
   IF( NB.LE.1 .OR. NB.GE.N ) THEN
!
!        Use unblocked code
!
      CALL DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
   ELSE
!
!        Use blocked code
!
      IF( UPPER ) THEN
!
!           Compute inverse of upper triangular matrix
!
         DO 20 J = 1, N, NB
            JB = MIN( NB, N-J+1 )
!
!              Compute rows 1:j-1 of current block column
!
            CALL DTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,&
                 JB, ONE, A, LDA, A( 1, J ), LDA )
            CALL DTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,&
                 JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
!
!              Compute inverse of current diagonal block
!
            CALL DTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
20       CONTINUE
      ELSE
!
!           Compute inverse of lower triangular matrix
!
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 30 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
            IF( J+JB.LE.N ) THEN
!
!                 Compute rows j+jb:n of current block column
!
               CALL DTRMM( 'Left', 'Lower', 'No transpose', DIAG,&
                    N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,&
                    A( J+JB, J ), LDA )
               CALL DTRSM( 'Right', 'Lower', 'No transpose', DIAG,&
                    N-J-JB+1, JB, -ONE, A( J, J ), LDA,&
                    A( J+JB, J ), LDA )
            END IF
!
!              Compute inverse of current diagonal block
!
            CALL DTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
30       CONTINUE
      END IF
   END IF
!
   RETURN
!
!     End of DTRTRI
!
 END SUBROUTINE DTRTRI
!
!=
!
!> \brief \b DTRTI2 computes the inverse of a triangular matrix (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DTRTI2 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrti2.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrti2.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrti2.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTRTI2 computes the inverse of a real upper or lower triangular
!> matrix.
!>
!> This is the Level 2 BLAS version of the algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the triangular matrix A.  If UPLO = 'U', the
!>          leading n by n upper triangular part of the array A contains
!>          the upper triangular matrix, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n by n lower triangular part of the array A contains
!>          the lower triangular matrix, and the strictly upper
!>          triangular part of A is not referenced.  If DIAG = 'U', the
!>          diagonal elements of A are also not referenced and are
!>          assumed to be 1.
!>
!>          On exit, the (triangular) inverse of the original matrix, in
!>          the same storage format.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
 SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
   CHARACTER          DIAG, UPLO
   INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE
   PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            NOUNIT, UPPER
   INTEGER            J
   DOUBLE PRECISION   AJJ
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DSCAL, DTRMV, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   NOUNIT = LSAME( DIAG, 'N' )
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
      INFO = -2
   ELSE IF( N.LT.0 ) THEN
      INFO = -3
   ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
      INFO = -5
   END IF
   IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'DTRTI2', -INFO )
      RETURN
   END IF
!
   IF( UPPER ) THEN
!
!        Compute inverse of upper triangular matrix.
!
      DO 10 J = 1, N
         IF( NOUNIT ) THEN
            A( J, J ) = ONE / A( J, J )
            AJJ = -A( J, J )
         ELSE
            AJJ = -ONE
         END IF
!
!           Compute elements 1:j-1 of j-th column.
!
         CALL DTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,&
              A( 1, J ), 1 )
         CALL DSCAL( J-1, AJJ, A( 1, J ), 1 )
10    CONTINUE
   ELSE
!
!        Compute inverse of lower triangular matrix.
!
      DO 20 J = N, 1, -1
         IF( NOUNIT ) THEN
            A( J, J ) = ONE / A( J, J )
            AJJ = -A( J, J )
         ELSE
            AJJ = -ONE
         END IF
         IF( J.LT.N ) THEN
!
!              Compute elements j+1:n of j-th column.
!
            CALL DTRMV( 'Lower', 'No transpose', DIAG, N-J,&
                 A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
            CALL DSCAL( N-J, AJJ, A( J+1, J ), 1 )
         END IF
20    CONTINUE
   END IF
!
   RETURN
!
!     End of DTRTI2
!
 END SUBROUTINE DTRTI2
!
!=
!
!> \brief \b DCOPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
! 
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DCOPY copies a vector, x, to a vector, y.
!>    uses unrolled loops for increments equal to one.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
 SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
    DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
    INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
    INTRINSIC MOD
!     ..
    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
       M = MOD(N,7)
       IF (M.NE.0) THEN
          DO I = 1,M
             DY(I) = DX(I)
          END DO
          IF (N.LT.7) RETURN
       END IF
       MP1 = M + 1
       DO I = MP1,N,7
          DY(I) = DX(I)
          DY(I+1) = DX(I+1)
          DY(I+2) = DX(I+2)
          DY(I+3) = DX(I+3)
          DY(I+4) = DX(I+4)
          DY(I+5) = DX(I+5)
          DY(I+6) = DX(I+6)
       END DO
    ELSE      
!
!        code for unequal increments or equal increments
!          not equal to 1
!
       IX = 1
       IY = 1
       IF (INCX.LT.0) IX = (-N+1)*INCX + 1
       IF (INCY.LT.0) IY = (-N+1)*INCY + 1
       DO I = 1,N
          DY(IY) = DX(IX)
          IX = IX + INCX
          IY = IY + INCY
       END DO
    END IF
    RETURN
  END SUBROUTINE DCOPY
!> \brief \b DDOT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
! 
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DDOT forms the dot product of two vectors.
!>    uses unrolled loops for increments equal to one.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DTEMP + DX(I)*DY(I)
            END DO
            IF (N.LT.5) THEN
               DDOT=DTEMP
            RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) + &
                 DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
         END DO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DTEMP + DX(IX)*DY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      DDOT = DTEMP
      RETURN
      END
!> \brief \b DGEMM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER K,LDA,LDB,LDC,M,N
!       CHARACTER TRANSA,TRANSB
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEMM  performs one of the matrix-matrix operations
!>
!>    C := alpha*op( A )*op( B ) + beta*C,
!>
!> where  op( X ) is one of
!>
!>    op( X ) = X   or   op( X ) = X**T,
!>
!> alpha and beta are scalars, and A, B and C are matrices, with op( A )
!> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>           On entry, TRANSA specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA = 'N' or 'n',  op( A ) = A.
!>
!>              TRANSA = 'T' or 't',  op( A ) = A**T.
!>
!>              TRANSA = 'C' or 'c',  op( A ) = A**T.
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER*1
!>           On entry, TRANSB specifies the form of op( B ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSB = 'N' or 'n',  op( B ) = B.
!>
!>              TRANSB = 'T' or 't',  op( B ) = B**T.
!>
!>              TRANSB = 'C' or 'c',  op( B ) = B**T.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry,  M  specifies  the number  of rows  of the  matrix
!>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry,  N  specifies the number  of columns of the matrix
!>           op( B ) and the number of columns of the matrix C. N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry,  K  specifies  the number of columns of the matrix
!>           op( A ) and the number of rows of the matrix op( B ). K must
!>           be at least  zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!>           part of the array  A  must contain the matrix  A,  otherwise
!>           the leading  k by m  part of the array  A  must contain  the
!>           matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!>           least  max( 1, k ).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!>           part of the array  B  must contain the matrix  B,  otherwise
!>           the leading  n by k  part of the array  B  must contain  the
!>           matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!>           least  max( 1, n ).
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION.
!>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!>           supplied as zero then C need not be set on input.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!>           Before entry, the leading  m by n  part of the array  C must
!>           contain the matrix  C,  except when  beta  is zero, in which
!>           case C need not be set on entry.
!>           On exit, the array  C  is overwritten by the  m by n  matrix
!>           ( alpha*op( A )*op( B ) + beta*C ).
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>           On entry, LDC specifies the first dimension of C as declared
!>           in  the  calling  (sub)  program.   LDC  must  be  at  least
!>           max( 1, m ).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2015
!
!> \ingroup double_blas_level3
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 3 Blas routine.
!>
!>  -- Written on 8-February-1989.
!>     Jack Dongarra, Argonne National Laboratory.
!>     Iain Duff, AERE Harwell.
!>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!>     Sven Hammarling, Numerical Algorithms Group Ltd.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!  -- Reference BLAS level3 routine (version 3.6.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!     ..
!
!  =====================================================================
!
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL NOTA,NOTB
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
          NROWA = M
          NCOLA = K
      ELSE
          NROWA = K
          NCOLA = M
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. &
           (.NOT.LSAME(TRANSA,'T'))) THEN
         INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND. &
           (.NOT.LSAME(TRANSB,'T'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMM ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
           (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
!
!     And if  alpha.eq.zero.
!
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
!
!     Start the operations.
!
      IF (NOTB) THEN
          IF (NOTA) THEN
!
!           Form  C := alpha*A*B + beta*C.
!
              DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                      TEMP = ALPHA*B(L,J)
                      DO 70 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE
!
!           Form  C := alpha*A**T*B + beta*C
!
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (NOTA) THEN
!
!           Form  C := alpha*A*B**T + beta*C
!
              DO 170 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 130 I = 1,M
                          C(I,J) = ZERO
  130                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 140 I = 1,M
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  END IF
                  DO 160 L = 1,K
                      TEMP = ALPHA*B(J,L)
                      DO 150 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  150                 CONTINUE
  160             CONTINUE
  170         CONTINUE
          ELSE
!
!           Form  C := alpha*A**T*B**T + beta*C
!
              DO 200 J = 1,N
                  DO 190 I = 1,M
                      TEMP = ZERO
                      DO 180 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
!
      RETURN
!
!     End of DGEMM .
!
      END
!> \brief \b DGEMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER INCX,INCY,LDA,M,N
!       CHARACTER TRANS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEMV  performs one of the matrix-vector operations
!>
!>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are vectors and A is an
!> m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!>
!>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!>
!>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of the matrix A.
!>           M must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, m ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array of DIMENSION at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!>           Before entry, the incremented array X must contain the
!>           vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION.
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array of DIMENSION at least
!>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!>           Before entry with BETA non-zero, the incremented array Y
!>           must contain the vector y. On exit, Y is overwritten by the
!>           updated vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2015
!
!> \ingroup double_blas_level2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!  -- Reference BLAS level2 routine (version 3.6.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
!     ..
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
           .NOT.LSAME(TRANS,'C')) THEN
         INFO = 1
      ELSE IF (M.LT.0) THEN
         INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMV ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
           ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
!
!        Form  y := alpha*A*x + y.
!
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  TEMP = ALPHA*X(JX)
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
!
!        Form  y := alpha*A**T*x + y.
!
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 100 J = 1,N
                  TEMP = ZERO
                  DO 90 I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  100         CONTINUE
          ELSE
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO 110 I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
!
      RETURN
!
!     End of DGEMV .
!
      END
!> \brief \b DGETRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DGETRF + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetrf.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetrf.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetrf.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGETRF computes an LU factorization of a general M-by-N matrix A
!> using partial pivoting with row interchanges.
!>
!> The factorization has the form
!>    A = P * L * U
!> where P is a permutation matrix, L is lower triangular with unit
!> diagonal elements (lower trapezoidal if m > n), and U is upper
!> triangular (upper trapezoidal if m < n).
!>
!> This is the right-looking Level 3 BLAS version of the algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the M-by-N matrix to be factored.
!>          On exit, the factors L and U from the factorization
!>          A = P*L*U; the unit diagonal elements of L are not stored.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (min(M,N))
!>          The pivot indices; for 1 <= i <= min(M,N), row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!>                has been completed, but the factor U is exactly
!>                singular, and division by zero will occur if it is used
!>                to solve a system of equations.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2015
!
!> \ingroup doubleGEcomputational
!
!  =====================================================================
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK computational routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEMM, DGETRF2, DLASWP, DTRSM, XERBLA
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL DGETRF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL DGETRF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 ) &
                 INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, &
                    IPIV, 1 )
!
!              Compute block row of U.
!
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
                    N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
                    LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                       N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, &
                       A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), &
                       LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of DGETRF
!
      END
!> \brief \b DGETRF2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       RECURSIVE SUBROUTINE DGETRF2( M, N, A, LDA, IPIV, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGETRF2 computes an LU factorization of a general M-by-N matrix A
!> using partial pivoting with row interchanges.
!>
!> The factorization has the form
!>    A = P * L * U
!> where P is a permutation matrix, L is lower triangular with unit
!> diagonal elements (lower trapezoidal if m > n), and U is upper
!> triangular (upper trapezoidal if m < n).
!>
!> This is the recursive version of the algorithm. It divides
!> the matrix into four submatrices:
!>            
!>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
!>    A = [ -----|----- ]  with n1 = min(m,n)
!>        [  A21 | A22  ]       n2 = n-n1
!>            
!>                                       [ A11 ]
!> The subroutine calls itself to factor [ --- ],
!>                                       [ A12 ]
!>                 [ A12 ]
!> do the swaps on [ --- ], solve A12, update A22,
!>                 [ A22 ]
!>
!> then calls itself to factor A22 and do the swaps on A21.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the M-by-N matrix to be factored.
!>          On exit, the factors L and U from the factorization
!>          A = P*L*U; the unit diagonal elements of L are not stored.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (min(M,N))
!>          The pivot indices; for 1 <= i <= min(M,N), row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!>                has been completed, but the factor U is exactly
!>                singular, and division by zero will occur if it is used
!>                to solve a system of equations.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2015
!
!> \ingroup doubleGEcomputational
!
!  =====================================================================
      RECURSIVE SUBROUTINE DGETRF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK computational routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN, TEMP
      INTEGER            I, IINFO, N1, N2
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      INTEGER            IDAMAX
!      EXTERNAL           DLAMCH, IDAMAX
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEMM, DSCAL, DLASWP, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      IF ( M.EQ.1 ) THEN
!
!        Use unblocked code for one row case
!        Just need to handle IPIV and INFO
!
         IPIV( 1 ) = 1
         IF ( A(1,1).EQ.ZERO ) INFO = 1
!
      ELSE IF( N.EQ.1 ) THEN
!
!        Use unblocked code for one column case
!
!
!        Compute machine safe minimum
!
         SFMIN = DLAMCH('S')
!
!        Find pivot and test for singularity
!
         I = IDAMAX( M, A( 1, 1 ), 1 )
         IPIV( 1 ) = I
         IF( A( I, 1 ).NE.ZERO ) THEN
!
!           Apply the interchange
!
            IF( I.NE.1 ) THEN
               TEMP = A( 1, 1 )
               A( 1, 1 ) = A( I, 1 )
               A( I, 1 ) = TEMP
            END IF
!
!           Compute elements 2:M of the column
!
            IF( ABS(A( 1, 1 )) .GE. SFMIN ) THEN
               CALL DSCAL( M-1, ONE / A( 1, 1 ), A( 2, 1 ), 1 )
            ELSE
               DO 10 I = 1, M-1
                  A( 1+I, 1 ) = A( 1+I, 1 ) / A( 1, 1 )
   10          CONTINUE
            END IF
!
         ELSE
            INFO = 1
         END IF
!
      ELSE
!
!        Use recursive code
!
         N1 = MIN( M, N ) / 2
         N2 = N-N1
!
!               [ A11 ]
!        Factor [ --- ]
!               [ A21 ]
!
         CALL DGETRF2( M, N1, A, LDA, IPIV, IINFO )

         IF ( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO
!
!                              [ A12 ]
!        Apply interchanges to [ --- ]
!                              [ A22 ]
!
         CALL DLASWP( N2, A( 1, N1+1 ), LDA, 1, N1, IPIV, 1 )
!
!        Solve A12
!
         CALL DTRSM( 'L', 'L', 'N', 'U', N1, N2, ONE, A, LDA,  &
              A( 1, N1+1 ), LDA )
!
!        Update A22
!
         CALL DGEMM( 'N', 'N', M-N1, N2, N1, -ONE, A( N1+1, 1 ), LDA,  &
                    A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA )
!
!        Factor A22
!
         CALL DGETRF2( M-N1, N2, A( N1+1, N1+1 ), LDA, IPIV( N1+1 ), &
              IINFO )
!
!        Adjust INFO and the pivot indices
!
         IF ( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + N1
         DO 20 I = N1+1, MIN( M, N )
            IPIV( I ) = IPIV( I ) + N1
   20    CONTINUE
!
!        Apply interchanges to A21
!
         CALL DLASWP( N1, A( 1, 1 ), LDA, N1+1, MIN( M, N), IPIV, 1 )
!
      END IF
      RETURN
!
!     End of DGETRF2
!
      END
!> \brief \b DGETRS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DGETRS + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetrs.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetrs.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetrs.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGETRS solves a system of linear equations
!>    A * X = B  or  A**T * X = B
!> with a general N-by-N matrix A using the LU factorization computed
!> by DGETRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations:
!>          = 'N':  A * X = B  (No transpose)
!>          = 'T':  A**T* X = B  (Transpose)
!>          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The factors L and U from the factorization A = P*L*U
!>          as computed by DGETRF.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!>          On entry, the right hand side matrix B.
!>          On exit, the solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup doubleGEcomputational
!
!  =====================================================================
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLASWP, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
           LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, &
              ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, &
              NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A**T * X = B.
!
!        Solve U**T *X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
              ONE, A, LDA, B, LDB )
!
!        Solve L**T *X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, &
              A, LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of DGETRS
!
    END SUBROUTINE DGETRS
!> \brief \b DISNAN tests input for NaN.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DISNAN + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/disnan.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/disnan.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/disnan.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION DISNAN( DIN )
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   DIN
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!> future.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIN
!> \verbatim
!>          DIN is DOUBLE PRECISION
!>          Input to test for NaN.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      LOGICAL FUNCTION DISNAN( DIN )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   DIN
!     ..
!
!=
!
!  .. External Functions ..
!      LOGICAL DLAISNAN
!      EXTERNAL DLAISNAN
!  ..
!  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)
      RETURN
      END
!> \brief \b DLAISNAN tests input for NaN by comparing two arguments for inequality.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAISNAN + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaisnan.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaisnan.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaisnan.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   DIN1, DIN2
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is not for general use.  It exists solely to avoid
!> over-optimization in DISNAN.
!>
!> DLAISNAN checks for NaNs by comparing its two arguments for
!> inequality.  NaN is the only floating-point value where NaN != NaN
!> returns .TRUE.  To check for NaNs, pass the same variable as both
!> arguments.
!>
!> A compiler must assume that the two arguments are
!> not the same variable, and the test will not be optimized away.
!> Interprocedural or whole-program optimization may delete this
!> test.  The ISNAN functions will be replaced by the correct
!> Fortran 03 intrinsic once the intrinsic is widely available.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIN1
!> \verbatim
!>          DIN1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] DIN2
!> \verbatim
!>          DIN2 is DOUBLE PRECISION
!>          Two numbers to compare for inequality.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   DIN1, DIN2
!     ..
!
!  =====================================================================
!
!  .. Executable Statements ..
      DLAISNAN = (DIN1.NE.DIN2)
      RETURN
      END
!> \brief \b DLAMCH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAMCH determines double precision machine parameters.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CMACH
!> \verbatim
!>          Specifies the value to be returned by DLAMCH:
!>          = 'E' or 'e',   DLAMCH := eps
!>          = 'S' or 's ,   DLAMCH := sfmin
!>          = 'B' or 'b',   DLAMCH := base
!>          = 'P' or 'p',   DLAMCH := eps*base
!>          = 'N' or 'n',   DLAMCH := t
!>          = 'R' or 'r',   DLAMCH := rnd
!>          = 'M' or 'm',   DLAMCH := emin
!>          = 'U' or 'u',   DLAMCH := rmin
!>          = 'L' or 'l',   DLAMCH := emax
!>          = 'O' or 'o',   DLAMCH := rmax
!>          where
!>          eps   = relative machine precision
!>          sfmin = safe minimum, such that 1/sfmin does not overflow
!>          base  = base of the machine
!>          prec  = eps*base
!>          t     = number of (base) digits in the mantissa
!>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!>          emin  = minimum exponent before (gradual) underflow
!>          rmin  = underflow threshold - base**(emin-1)
!>          emax  = largest exponent before overflow
!>          rmax  = overflow threshold  - (base**emax)*(1-eps)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2015
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!
!  -- LAPACK auxiliary routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   RND, EPS, SFMIN, SMALL, RMACH
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT, &
           MINEXPONENT, RADIX, TINY
!     ..
!     .. Executable Statements ..
!
!
!     Assume rounding, not chopping. Always.
!
      RND = ONE
!
      IF( ONE.EQ.RND ) THEN
         EPS = EPSILON(ZERO) * 0.5
      ELSE
         EPS = EPSILON(ZERO)
      END IF
!
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         SFMIN = TINY(ZERO)
         SMALL = ONE / HUGE(ZERO)
         IF( SMALL.GE.SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            SFMIN = SMALL*( ONE+EPS )
         END IF
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = EPS * RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = DIGITS(ZERO)
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = MINEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = tiny(zero)
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = MAXEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = HUGE(ZERO)
      ELSE
         RMACH = ZERO
      END IF
!
      DLAMCH = RMACH
      RETURN
!
!     End of DLAMCH
!
      END
!***********************************************************************
!> \brief \b DLAMC3
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!> the addition of  A  and  B ,  for use in situations where optimizers
!> might hold one of these in a register.
!> \endverbatim
!> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!> \date November 2015
!> \ingroup auxOTHERauxiliary
!>
!> \param[in] A
!> \verbatim
!>          A is a DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is a DOUBLE PRECISION
!>          The values A and B.
!> \endverbatim
!>
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
!
!  -- LAPACK auxiliary routine (version 3.6.0) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2010
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
!     ..
! =====================================================================
!
!     .. Executable Statements ..
!
      DLAMC3 = A + B
!
      RETURN
!
!     End of DLAMC3
!
      END
!
!***********************************************************************
!> \brief \b DLAMRG creates a permutation list to merge the entries of two independently sorted sets into a single set sorted in ascending order.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAMRG + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlamrg.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlamrg.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlamrg.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAMRG( N1, N2, A, DTRD1, DTRD2, INDEX )
! 
!       .. Scalar Arguments ..
!       INTEGER            DTRD1, DTRD2, N1, N2
!       ..
!       .. Array Arguments ..
!       INTEGER            INDEX( * )
!       DOUBLE PRECISION   A( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAMRG will create a permutation list which will merge the elements
!> of A (which is composed of two independently sorted sets) into a
!> single set which is sorted in ascending order.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!>         These arguements contain the respective lengths of the two
!>         sorted lists to be merged.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (N1+N2)
!>         The first N1 elements of A contain a list of numbers which
!>         are sorted in either ascending or descending order.  Likewise
!>         for the final N2 elements.
!> \endverbatim
!>
!> \param[in] DTRD1
!> \verbatim
!>          DTRD1 is INTEGER
!> \endverbatim
!>
!> \param[in] DTRD2
!> \verbatim
!>          DTRD2 is INTEGER
!>         These are the strides to be taken through the array A.
!>         Allowable strides are 1 and -1.  They indicate whether a
!>         subset of A is sorted in ascending (DTRDx = 1) or descending
!>         (DTRDx = -1) order.
!> \endverbatim
!>
!> \param[out] INDEX
!> \verbatim
!>          INDEX is INTEGER array, dimension (N1+N2)
!>         On exit this array will contain a permutation such that
!>         if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be
!>         sorted in ascending order.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DLAMRG( N1, N2, A, DTRD1, DTRD2, INDEX )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            DTRD1, DTRD2, N1, N2
!     ..
!     .. Array Arguments ..
      INTEGER            INDEX( * )
      DOUBLE PRECISION   A( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IND1, IND2, N1SV, N2SV
!     ..
!     .. Executable Statements ..
!
      N1SV = N1
      N2SV = N2
      IF( DTRD1.GT.0 ) THEN
         IND1 = 1
      ELSE
         IND1 = N1
      END IF
      IF( DTRD2.GT.0 ) THEN
         IND2 = 1 + N1
      ELSE
         IND2 = N1 + N2
      END IF
      I = 1
!     while ( (N1SV > 0) & (N2SV > 0) )
   10 CONTINUE
      IF( N1SV.GT.0 .AND. N2SV.GT.0 ) THEN
         IF( A( IND1 ).LE.A( IND2 ) ) THEN
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + DTRD1
            N1SV = N1SV - 1
         ELSE
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + DTRD2
            N2SV = N2SV - 1
         END IF
         GO TO 10
      END IF
!     end while
      IF( N1SV.EQ.0 ) THEN
         DO 20 N1SV = 1, N2SV
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + DTRD2
   20    CONTINUE
      ELSE
!     N2SV .EQ. 0
         DO 30 N2SV = 1, N1SV
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + DTRD1
   30    CONTINUE
      END IF
!
      RETURN
!
!     End of DLAMRG
!
      END
!> \brief \b DLASWP performs a series of row interchanges on a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLASWP + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaswp.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaswp.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaswp.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
! 
!       .. Scalar Arguments ..
!       INTEGER            INCX, K1, K2, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASWP performs a series of row interchanges on the matrix A.
!> One row interchange is initiated for each of rows K1 through K2 of A.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the matrix of column dimension N to which the row
!>          interchanges will be applied.
!>          On exit, the permuted matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!> \endverbatim
!>
!> \param[in] K1
!> \verbatim
!>          K1 is INTEGER
!>          The first element of IPIV for which a row interchange will
!>          be done.
!> \endverbatim
!>
!> \param[in] K2
!> \verbatim
!>          K2 is INTEGER
!>          The last element of IPIV for which a row interchange will
!>          be done.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (K2*abs(INCX))
!>          The vector of pivot indices.  Only the elements in positions
!>          K1 through K2 of IPIV are accessed.
!>          IPIV(K) = L implies rows K and L are to be interchanged.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of IPIV.  If IPIV
!>          is negative, the pivots are applied in reverse order.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup doubleOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Modified by
!>   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
!
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
!
      RETURN
!
!     End of DLASWP
!
      END
!> \brief \b DLASYF computes a partial factorization of a real symmetric matrix using the Bunch-Kaufman diagonal pivoting method.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASYF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasyf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasyf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasyf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KB, LDA, LDW, N, NB
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), W( LDW, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASYF computes a partial factorization of a real symmetric matrix A
!> using the Bunch-Kaufman diagonal pivoting method. The partial
!> factorization has the form:
!>
!> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
!>       ( 0  U22 ) (  0   D  ) ( U12**T U22**T )
!>
!> A  =  ( L11  0 ) (  D   0  ) ( L11**T L21**T )  if UPLO = 'L'
!>       ( L21  I ) (  0  A22 ) (  0       I    )
!>
!> where the order of D is at most NB. The actual order is returned in
!> the argument KB, and is either NB or NB-1, or N if N <= NB.
!>
!> DLASYF is an auxiliary routine called by DSYTRF. It uses blocked code
!> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
!> A22 (if UPLO = 'L').
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is stored:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The maximum number of columns of the matrix A that should be
!>          factored.  NB should be at least 2 to allow for 2-by-2 pivot
!>          blocks.
!> \endverbatim
!>
!> \param[out] KB
!> \verbatim
!>          KB is INTEGER
!>          The number of columns of A that were actually factored.
!>          KB is either NB-1 or NB, or N if N <= NB.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          n-by-n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n-by-n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>          On exit, A contains details of the partial factorization.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D.
!>
!>          If UPLO = 'U':
!>             Only the last KB elements of IPIV are set.
!>
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>             interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) = IPIV(k-1) < 0, then rows and columns
!>             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!>             is a 2-by-2 diagonal block.
!>
!>          If UPLO = 'L':
!>             Only the first KB elements of IPIV are set.
!>
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>             interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) = IPIV(k+1) < 0, then rows and columns
!>             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1)
!>             is a 2-by-2 diagonal block.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (LDW,NB)
!> \endverbatim
!>
!> \param[in] LDW
!> \verbatim
!>          LDW is INTEGER
!>          The leading dimension of the array W.  LDW >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!>               has been completed, but the block diagonal matrix D is
!>               exactly singular.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2013
!
!> \ingroup doubleSYcomputational
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  November 2013,  Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE DLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
!
!  -- LAPACK computational routine (version 3.5.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2013
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KB, LDA, LDW, N, NB
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), W( LDW, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP, &
           KSTEP, KW
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D21, D22, R1, &
           ROWMAX, T
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            IDAMAX
!      EXTERNAL           LSAME, IDAMAX
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DGEMM, DGEMV, DSCAL, DSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Initialize ALPHA for use in choosing pivot block size.
!
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Factorize the trailing columns of A using the upper triangle
!        of A and working backwards, and compute the matrix W = U12*D
!        for use in updating A11
!
!        K is the main loop index, decreasing from N in steps of 1 or 2
!
!        KW is the column of W which corresponds to column K of A
!
         K = N
   10    CONTINUE
         KW = NB + K - N
!
!        Exit from loop
!
         IF( ( K.LE.N-NB+1 .AND. NB.LT.N ) .OR. K.LT.1 ) GO TO 30
!
!        Copy column K of A to column KW of W and update it
!
         CALL DCOPY( K, A( 1, K ), 1, W( 1, KW ), 1 )
         IF( K.LT.N ) &
           CALL DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), LDA, &
                       W( K, KW+1 ), LDW, ONE, W( 1, KW ), 1 )
!
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = ABS( W( K, KW ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
         IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, W( 1, KW ), 1 )
            COLMAX = ABS( W( IMAX, KW ) )
         ELSE
            COLMAX = ZERO
         END IF
!
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
!
!           Column K is zero or underflow: set INFO and continue
!
            IF( INFO.EQ.0 ) INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            ELSE
!
!              Copy column IMAX to column KW-1 of W and update it
!
               CALL DCOPY( IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 )
               CALL DCOPY( K-IMAX, A( IMAX, IMAX+1 ), LDA, &
                    W( IMAX+1, KW-1 ), 1 )
               IF( K.LT.N ) &
                    CALL DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), &
                    LDA, W( IMAX, KW+1 ), LDW, ONE, &
                    W( 1, KW-1 ), 1 )
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = IMAX + IDAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 )
               ROWMAX = ABS( W( JMAX, KW-1 ) )
               IF( IMAX.GT.1 ) THEN
                  JMAX = IDAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( W( JMAX, KW-1 ) ) )
               END IF
!
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               ELSE IF( ABS( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
!
!                 copy column KW-1 of W to column KW of W
!
                  CALL DCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
               ELSE
!
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
!
!           ============================================================
!
!           KK is the column of A where pivoting step stopped
!
            KK = K - KSTEP + 1
!
!           KKW is the column of W which corresponds to column KK of A
!
            KKW = NB + KK - N
!
!           Interchange rows and columns KP and KK.
!           Updated column KP is already stored in column KKW of W.
!
            IF( KP.NE.KK ) THEN
!
!              Copy non-updated column KK to column KP of submatrix A
!              at step K. No need to copy element into column K
!              (or K and K-1 for 2-by-2 pivot) of A, since these columns
!              will be later overwritten.
!
               A( KP, KP ) = A( KK, KK )
               CALL DCOPY( KK-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ), &
                    LDA )
               IF( KP.GT.1 ) &
                 CALL DCOPY( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
!
!              Interchange rows KK and KP in last K+1 to N columns of A
!              (columns K (or K and K-1 for 2-by-2 pivot) of A will be
!              later overwritten). Interchange rows KK and KP
!              in last KKW to NB columns of W.
!
               IF( K.LT.N ) &
                    CALL DSWAP( N-K, A( KK, K+1 ), LDA, A( KP, K+1 ), &
                    LDA )
               CALL DSWAP( N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ), &
                    LDW )
            END IF
!
            IF( KSTEP.EQ.1 ) THEN
!
!              1-by-1 pivot block D(k): column kw of W now holds
!
!              W(kw) = U(k)*D(k),
!
!              where U(k) is the k-th column of U
!
!              Store subdiag. elements of column U(k)
!              and 1-by-1 block D(k) in column k of A.
!              NOTE: Diagonal element U(k,k) is a UNIT element
!              and not stored.
!                 A(k,k) := D(k,k) = W(k,kw)
!                 A(1:k-1,k) := U(1:k-1,k) = W(1:k-1,kw)/D(k,k)
!
               CALL DCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
               R1 = ONE / A( K, K )
               CALL DSCAL( K-1, R1, A( 1, K ), 1 )
!
            ELSE
!
!              2-by-2 pivot block D(k): columns kw and kw-1 of W now hold
!
!              ( W(kw-1) W(kw) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
!              Store U(1:k-2,k-1) and U(1:k-2,k) and 2-by-2
!              block D(k-1:k,k-1:k) in columns k-1 and k of A.
!              NOTE: 2-by-2 diagonal block U(k-1:k,k-1:k) is a UNIT
!              block and not stored.
!                 A(k-1:k,k-1:k) := D(k-1:k,k-1:k) = W(k-1:k,kw-1:kw)
!                 A(1:k-2,k-1:k) := U(1:k-2,k:k-1:k) =
!                 = W(1:k-2,kw-1:kw) * ( D(k-1:k,k-1:k)**(-1) )
!
               IF( K.GT.2 ) THEN
!
!                 Compose the columns of the inverse of 2-by-2 pivot
!                 block D in the following way to reduce the number
!                 of FLOPS when we myltiply panel ( W(kw-1) W(kw) ) by
!                 this inverse
!
!                 D**(-1) = ( d11 d21 )**(-1) =
!                           ( d21 d22 )
!
!                 = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
!                                        ( (-d21 ) ( d11 ) )
!
!                 = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
!
!                   * ( ( d22/d21 ) (      -1 ) ) =
!                     ( (      -1 ) ( d11/d21 ) )
!
!                 = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) (  -1 ) ) =
!                                           ( ( -1  ) ( D22 ) )
!
!                 = 1/d21 * T * ( ( D11 ) (  -1 ) )
!                               ( (  -1 ) ( D22 ) )
!
!                 = D21 * ( ( D11 ) (  -1 ) )
!                         ( (  -1 ) ( D22 ) )
!
                  D21 = W( K-1, KW )
                  D11 = W( K, KW ) / D21
                  D22 = W( K-1, KW-1 ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
!
!                 Update elements in columns A(k-1) and A(k) as
!                 dot products of rows of ( W(kw-1) W(kw) ) and columns
!                 of D**(-1)
!
                  DO 20 J = 1, K - 2
                     A( J, K-1 ) = D21*( D11*W( J, KW-1 )-W( J, KW ) )
                     A( J, K ) = D21*( D22*W( J, KW )-W( J, KW-1 ) )
   20             CONTINUE
               END IF
!
!              Copy D(k) to A
!
               A( K-1, K-1 ) = W( K-1, KW-1 )
               A( K-1, K ) = W( K-1, KW )
               A( K, K ) = W( K, KW )
!
            END IF
!
         END IF
!
!        Store details of the interchanges in IPIV
!
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
!
!        Decrease K and return to the start of the main loop
!
         K = K - KSTEP
         GO TO 10
!
   30    CONTINUE
!
!        Update the upper triangle of A11 (= A(1:k,1:k)) as
!
!        A11 := A11 - U12*D*U12**T = A11 - U12*W**T
!
!        computing blocks of NB columns at a time
!
         DO 50 J = ( ( K-1 ) / NB )*NB + 1, 1, -NB
            JB = MIN( NB, K-J+1 )
!
!           Update the upper triangle of the diagonal block
!
            DO 40 JJ = J, J + JB - 1
               CALL DGEMV( 'No transpose', JJ-J+1, N-K, -ONE, &
                    A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, ONE, &
                    A( J, JJ ), 1 )
   40       CONTINUE
!
!           Update the rectangular superdiagonal block
!
            CALL DGEMM( 'No transpose', 'Transpose', J-1, JB, N-K, -ONE, &
                 A( 1, K+1 ), LDA, W( J, KW+1 ), LDW, ONE, &
                 A( 1, J ), LDA )
   50    CONTINUE
!
!        Put U12 in standard form by partially undoing the interchanges
!        in columns k+1:n looping backwards from k+1 to n
!
         J = K + 1
   60    CONTINUE
!
!           Undo the interchanges (if any) of rows JJ and JP at each
!           step J
!
!           (Here, J is a diagonal index)
            JJ = J
            JP = IPIV( J )
            IF( JP.LT.0 ) THEN
               JP = -JP
!              (Here, J is a diagonal index)
               J = J + 1
            END IF
!           (NOTE: Here, J is used to determine row length. Length N-J+1
!           of the rows to swap back doesn't include diagonal element)
            J = J + 1
            IF( JP.NE.JJ .AND. J.LE.N ) &
              CALL DSWAP( N-J+1, A( JP, J ), LDA, A( JJ, J ), LDA )
         IF( J.LT.N ) GO TO 60
!
!        Set KB to the number of columns factorized
!
         KB = N - K
!
      ELSE
!
!        Factorize the leading columns of A using the lower triangle
!        of A and working forwards, and compute the matrix W = L21*D
!        for use in updating A22
!
!        K is the main loop index, increasing from 1 in steps of 1 or 2
!
         K = 1
   70    CONTINUE
!
!        Exit from loop
!
         IF( ( K.GE.NB .AND. NB.LT.N ) .OR. K.GT.N ) GO TO 90
!
!        Copy column K of A to column K of W and update it
!
         CALL DCOPY( N-K+1, A( K, K ), 1, W( K, K ), 1 )
         CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ), LDA, &
              W( K, 1 ), LDW, ONE, W( K, K ), 1 )
!
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = ABS( W( K, K ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
         IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, W( K+1, K ), 1 )
            COLMAX = ABS( W( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
!
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
!
!           Column K is zero or underflow: set INFO and continue
!
            IF( INFO.EQ.0 ) INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            ELSE
!
!              Copy column IMAX to column K+1 of W and update it
!
               CALL DCOPY( IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1 )
               CALL DCOPY( N-IMAX+1, A( IMAX, IMAX ), 1, W( IMAX, K+1 ), &
                    1 )
               CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ), &
                    LDA, W( IMAX, 1 ), LDW, ONE, W( K, K+1 ), 1 )
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = K - 1 + IDAMAX( IMAX-K, W( K, K+1 ), 1 )
               ROWMAX = ABS( W( JMAX, K+1 ) )
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + IDAMAX( N-IMAX, W( IMAX+1, K+1 ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( W( JMAX, K+1 ) ) )
               END IF
!
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               ELSE IF( ABS( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
!
!                 copy column K+1 of W to column K of W
!
                  CALL DCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
               ELSE
!
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
!
!           ============================================================
!
!           KK is the column of A where pivoting step stopped
!
            KK = K + KSTEP - 1
!
!           Interchange rows and columns KP and KK.
!           Updated column KP is already stored in column KK of W.
!
            IF( KP.NE.KK ) THEN
!
!              Copy non-updated column KK to column KP of submatrix A
!              at step K. No need to copy element into column K
!              (or K and K+1 for 2-by-2 pivot) of A, since these columns
!              will be later overwritten.
!
               A( KP, KP ) = A( KK, KK )
               CALL DCOPY( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), &
                    LDA )
               IF( KP.LT.N ) &
                    CALL DCOPY( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
!
!              Interchange rows KK and KP in first K-1 columns of A
!              (columns K (or K and K+1 for 2-by-2 pivot) of A will be
!              later overwritten). Interchange rows KK and KP
!              in first KK columns of W.
!
               IF( K.GT.1 ) &
                    CALL DSWAP( K-1, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
               CALL DSWAP( KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW )
            END IF
!
            IF( KSTEP.EQ.1 ) THEN
!
!              1-by-1 pivot block D(k): column k of W now holds
!
!              W(k) = L(k)*D(k),
!
!              where L(k) is the k-th column of L
!
!              Store subdiag. elements of column L(k)
!              and 1-by-1 block D(k) in column k of A.
!              (NOTE: Diagonal element L(k,k) is a UNIT element
!              and not stored)
!                 A(k,k) := D(k,k) = W(k,k)
!                 A(k+1:N,k) := L(k+1:N,k) = W(k+1:N,k)/D(k,k)
!
               CALL DCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
               IF( K.LT.N ) THEN
                  R1 = ONE / A( K, K )
                  CALL DSCAL( N-K, R1, A( K+1, K ), 1 )
               END IF
!
            ELSE
!
!              2-by-2 pivot block D(k): columns k and k+1 of W now hold
!
!              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
!
!              where L(k) and L(k+1) are the k-th and (k+1)-th columns
!              of L
!
!              Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2
!              block D(k:k+1,k:k+1) in columns k and k+1 of A.
!              (NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT
!              block and not stored)
!                 A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1)
!                 A(k+2:N,k:k+1) := L(k+2:N,k:k+1) =
!                 = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) )
!
               IF( K.LT.N-1 ) THEN
!
!                 Compose the columns of the inverse of 2-by-2 pivot
!                 block D in the following way to reduce the number
!                 of FLOPS when we myltiply panel ( W(k) W(k+1) ) by
!                 this inverse
!
!                 D**(-1) = ( d11 d21 )**(-1) =
!                           ( d21 d22 )
!
!                 = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
!                                        ( (-d21 ) ( d11 ) )
!
!                 = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
!
!                   * ( ( d22/d21 ) (      -1 ) ) =
!                     ( (      -1 ) ( d11/d21 ) )
!
!                 = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) (  -1 ) ) =
!                                           ( ( -1  ) ( D22 ) )
!
!                 = 1/d21 * T * ( ( D11 ) (  -1 ) )
!                               ( (  -1 ) ( D22 ) )
!
!                 = D21 * ( ( D11 ) (  -1 ) )
!                         ( (  -1 ) ( D22 ) )
!
                  D21 = W( K+1, K )
                  D11 = W( K+1, K+1 ) / D21
                  D22 = W( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
!
!                 Update elements in columns A(k) and A(k+1) as
!                 dot products of rows of ( W(k) W(k+1) ) and columns
!                 of D**(-1)
!
                  DO 80 J = K + 2, N
                     A( J, K ) = D21*( D11*W( J, K )-W( J, K+1 ) )
                     A( J, K+1 ) = D21*( D22*W( J, K+1 )-W( J, K ) )
   80             CONTINUE
               END IF
!
!              Copy D(k) to A
!
               A( K, K ) = W( K, K )
               A( K+1, K ) = W( K+1, K )
               A( K+1, K+1 ) = W( K+1, K+1 )
!
            END IF
!
         END IF
!
!        Store details of the interchanges in IPIV
!
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
!
!        Increase K and return to the start of the main loop
!
         K = K + KSTEP
         GO TO 70
!
   90    CONTINUE
!
!        Update the lower triangle of A22 (= A(k:n,k:n)) as
!
!        A22 := A22 - L21*D*L21**T = A22 - L21*W**T
!
!        computing blocks of NB columns at a time
!
         DO 110 J = K, N, NB
            JB = MIN( NB, N-J+1 )
!
!           Update the lower triangle of the diagonal block
!
            DO 100 JJ = J, J + JB - 1
               CALL DGEMV( 'No transpose', J+JB-JJ, K-1, -ONE, &
                    A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, ONE, &
                    A( JJ, JJ ), 1 )
  100       CONTINUE
!
!           Update the rectangular subdiagonal block
!
            IF( J+JB.LE.N ) &
                 CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB, &
                 K-1, -ONE, A( J+JB, 1 ), LDA, W( J, 1 ), LDW, &
                 ONE, A( J+JB, J ), LDA )
  110    CONTINUE
!
!        Put L21 in standard form by partially undoing the interchanges
!        of rows in columns 1:k-1 looping backwards from k-1 to 1
!
         J = K - 1
  120    CONTINUE
!
!           Undo the interchanges (if any) of rows JJ and JP at each
!           step J
!
!           (Here, J is a diagonal index)
            JJ = J
            JP = IPIV( J )
            IF( JP.LT.0 ) THEN
               JP = -JP
!              (Here, J is a diagonal index)
               J = J - 1
            END IF
!           (NOTE: Here, J is used to determine row length. Length J
!           of the rows to swap back doesn't include diagonal element)
            J = J - 1
            IF( JP.NE.JJ .AND. J.GE.1 ) &
                 CALL DSWAP( J, A( JP, 1 ), LDA, A( JJ, 1 ), LDA )
         IF( J.GT.1 ) GO TO 120
!
!        Set KB to the number of columns factorized
!
         KB = K - 1
!
      END IF
      RETURN
!
!     End of DLASYF
!
      END
!> \brief \b DSCAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSCAL(N,DA,DX,INCX)
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION DA
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DSCAL scales a vector by a constant.
!>    uses unrolled loops for increment equal to one.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DSCAL(N,DA,DX,INCX)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DX(I) = DA*DX(I)
            END DO
            IF (N.LT.5) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            DX(I) = DA*DX(I)
            DX(I+1) = DA*DX(I+1)
            DX(I+2) = DA*DX(I+2)
            DX(I+3) = DA*DX(I+3)
            DX(I+4) = DA*DX(I+4)
         END DO
      ELSE
!
!        code for increment not equal to 1
!
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DX(I) = DA*DX(I)
         END DO
      END IF
      RETURN
      END
!> \brief \b DSWAP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
! 
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    interchanges two vectors.
!>    uses unrolled loops for increments equal one.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
         M = MOD(N,3)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DX(I)
               DX(I) = DY(I)
               DY(I) = DTEMP
            END DO
            IF (N.LT.3) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,3
            DTEMP = DX(I)
            DX(I) = DY(I)
            DY(I) = DTEMP
            DTEMP = DX(I+1)
            DX(I+1) = DY(I+1)
            DY(I+1) = DTEMP
            DTEMP = DX(I+2)
            DX(I+2) = DY(I+2)
            DY(I+2) = DTEMP
         END DO
      ELSE
!
!       code for unequal increments or equal increments not equal
!         to 1
!
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DX(IX)
            DX(IX) = DY(IY)
            DY(IY) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
!> \brief \b DSYMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER INCX,INCY,LDA,N
!       CHARACTER UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYMV  performs the matrix-vector  operation
!>
!>    y := alpha*A*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are n element vectors and
!> A is an n by n symmetric matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the array A is to be referenced as
!>           follows:
!>
!>              UPLO = 'U' or 'u'   Only the upper triangular part of A
!>                                  is to be referenced.
!>
!>              UPLO = 'L' or 'l'   Only the lower triangular part of A
!>                                  is to be referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!>           Before entry with  UPLO = 'U' or 'u', the leading n by n
!>           upper triangular part of the array A must contain the upper
!>           triangular part of the symmetric matrix and the strictly
!>           lower triangular part of A is not referenced.
!>           Before entry with UPLO = 'L' or 'l', the leading n by n
!>           lower triangular part of the array A must contain the lower
!>           triangular part of the symmetric matrix and the strictly
!>           upper triangular part of A is not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, n ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array of dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
!>           element vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION.
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array of dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y. On exit, Y is overwritten by the updated
!>           vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!  -- Reference BLAS level2 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
!     ..
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 5
      ELSE IF (INCX.EQ.0) THEN
          INFO = 7
      ELSE IF (INCY.EQ.0) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYMV ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!
!     Set up the start points in  X  and  Y.
!
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (N-1)*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
!     First form  y := beta*y.
!
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,N
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,N
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,N
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,N
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(UPLO,'U')) THEN
!
!        Form  y  when A is stored in upper triangle.
!
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  DO 50 I = 1,J - 1
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   50             CONTINUE
                  Y(J) = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2
   60         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 80 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  IX = KX
                  IY = KY
                  DO 70 I = 1,J - 1
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
                      IX = IX + INCX
                      IY = IY + INCY
   70             CONTINUE
                  Y(JY) = Y(JY) + TEMP1*A(J,J) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
   80         CONTINUE
          END IF
      ELSE
!
!        Form  y  when A is stored in lower triangle.
!
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 100 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  Y(J) = Y(J) + TEMP1*A(J,J)
                  DO 90 I = J + 1,N
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   90             CONTINUE
                  Y(J) = Y(J) + ALPHA*TEMP2
  100         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 120 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  Y(JY) = Y(JY) + TEMP1*A(J,J)
                  IX = JX
                  IY = JY
                  DO 110 I = J + 1,N
                      IX = IX + INCX
                      IY = IY + INCY
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
!
      RETURN
!
!     End of DSYMV .
!
      END
!> \brief \b DSYR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA
!       INTEGER INCX,LDA,N
!       CHARACTER UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),X(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYR   performs the symmetric rank 1 operation
!>
!>    A := alpha*x*x**T + A,
!>
!> where alpha is a real scalar, x is an n element vector and A is an
!> n by n symmetric matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the array A is to be referenced as
!>           follows:
!>
!>              UPLO = 'U' or 'u'   Only the upper triangular part of A
!>                                  is to be referenced.
!>
!>              UPLO = 'L' or 'l'   Only the lower triangular part of A
!>                                  is to be referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array of dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
!>           element vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!>           Before entry with  UPLO = 'U' or 'u', the leading n by n
!>           upper triangular part of the array A must contain the upper
!>           triangular part of the symmetric matrix and the strictly
!>           lower triangular part of A is not referenced. On exit, the
!>           upper triangular part of the array A is overwritten by the
!>           upper triangular part of the updated matrix.
!>           Before entry with UPLO = 'L' or 'l', the leading n by n
!>           lower triangular part of the array A must contain the lower
!>           triangular part of the symmetric matrix and the strictly
!>           upper triangular part of A is not referenced. On exit, the
!>           lower triangular part of the array A is overwritten by the
!>           lower triangular part of the updated matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, n ).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
!
!  -- Reference BLAS level2 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER INCX,LDA,N
      CHARACTER UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
!     ..
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 7
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYR  ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!
!     Set the start point in X if the increment is not unity.
!
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
      IF (LSAME(UPLO,'U')) THEN
!
!        Form  A  when A is stored in upper triangle.
!
          IF (INCX.EQ.1) THEN
              DO 20 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*X(J)
                      DO 10 I = 1,J
                          A(I,J) = A(I,J) + X(I)*TEMP
   10                 CONTINUE
                  END IF
   20         CONTINUE
          ELSE
              JX = KX
              DO 40 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = KX
                      DO 30 I = 1,J
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                  END IF
                  JX = JX + INCX
   40         CONTINUE
          END IF
      ELSE
!
!        Form  A  when A is stored in lower triangle.
!
          IF (INCX.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*X(J)
                      DO 50 I = J,N
                          A(I,J) = A(I,J) + X(I)*TEMP
   50                 CONTINUE
                  END IF
   60         CONTINUE
          ELSE
              JX = KX
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = JX
                      DO 70 I = J,N
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
   80         CONTINUE
          END IF
      END IF
!
      RETURN
!
!     End of DSYR  .
!
      END
!> \brief \b DSYTF2 computes the factorization of a real symmetric indefinite matrix, using the diagonal pivoting method (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DSYTF2 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytf2.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytf2.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytf2.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYTF2( UPLO, N, A, LDA, IPIV, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYTF2 computes the factorization of a real symmetric matrix A using
!> the Bunch-Kaufman diagonal pivoting method:
!>
!>    A = U*D*U**T  or  A = L*D*L**T
!>
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, U**T is the transpose of U, and D is symmetric and
!> block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> This is the unblocked version of the algorithm, calling Level 2 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is stored:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          n-by-n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n-by-n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, the block diagonal matrix D and the multipliers used
!>          to obtain the factor U or L (see below for further details).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D.
!>
!>          If UPLO = 'U':
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>             interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) = IPIV(k-1) < 0, then rows and columns
!>             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!>             is a 2-by-2 diagonal block.
!>
!>          If UPLO = 'L':
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>             interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) = IPIV(k+1) < 0, then rows and columns
!>             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1)
!>             is a 2-by-2 diagonal block.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
!>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!>               has been completed, but the block diagonal matrix D is
!>               exactly singular, and division by zero will occur if it
!>               is used to solve a system of equations.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2013
!
!> \ingroup doubleSYcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', then A = U*D*U**T, where
!>     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!>
!>             (   I    v    0   )   k-s
!>     U(k) =  (   0    I    0   )   s
!>             (   0    0    I   )   n-k
!>                k-s   s   n-k
!>
!>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!>  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!>
!>  If UPLO = 'L', then A = L*D*L**T, where
!>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!>
!>             (   I    0     0   )  k-1
!>     L(k) =  (   0    I     0   )  s
!>             (   0    v     I   )  n-k-s+1
!>                k-1   s  n-k-s+1
!>
!>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  09-29-06 - patch from
!>    Bobby Cheng, MathWorks
!>
!>    Replace l.204 and l.372
!>         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
!>    by
!>         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN
!>
!>  01-01-96 - Based on modifications by
!>    J. Lewis, Boeing Computer Services Company
!>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!>  1-96 - Based on modifications by J. Lewis, Boeing Computer Services
!>         Company
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE DSYTF2( UPLO, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK computational routine (version 3.5.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2013
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IMAX, J, JMAX, K, KK, KP, KSTEP
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22, R1, &
                        ROWMAX, T, WK, WKM1, WKP1
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME, DISNAN
!      INTEGER            IDAMAX
!      EXTERNAL           LSAME, IDAMAX, DISNAN
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DSCAL, DSWAP, DSYR, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTF2', -INFO )
         RETURN
      END IF
!
!     Initialize ALPHA for use in choosing pivot block size.
!
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
!
      IF( UPPER ) THEN
!
!        Factorize A as U*D*U**T using the upper triangle of A
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2
!
         K = N
   10    CONTINUE
!
!        If K < 1, exit from loop
!
         IF( K.LT.1 ) GO TO 70
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = ABS( A( K, K ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
         IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, A( 1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
!
         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN
!
!           Column K is zero or underflow, or contains a NaN:
!           set INFO and continue
!
            IF( INFO.EQ.0 ) INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            ELSE
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = IMAX + IDAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
               ROWMAX = ABS( A( IMAX, JMAX ) )
               IF( IMAX.GT.1 ) THEN
                  JMAX = IDAMAX( IMAX-1, A( 1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( A( JMAX, IMAX ) ) )
               END IF
!
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               ELSE IF( ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
               ELSE
!
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
!
            KK = K - KSTEP + 1
            IF( KP.NE.KK ) THEN
!
!              Interchange rows and columns KK and KP in the leading
!              submatrix A(1:k,1:k)
!
               CALL DSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
               CALL DSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ), &
                    LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
            END IF
!
!           Update the leading submatrix
!
            IF( KSTEP.EQ.1 ) THEN
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = U(k)*D(k)
!
!              where U(k) is the k-th column of U
!
!              Perform a rank-1 update of A(1:k-1,1:k-1) as
!
!              A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T
!
               R1 = ONE / A( K, K )
               CALL DSYR( UPLO, K-1, -R1, A( 1, K ), 1, A, LDA )
!
!              Store U(k) in column k
!
               CALL DSCAL( K-1, R1, A( 1, K ), 1 )
            ELSE
!
!              2-by-2 pivot block D(k): columns k and k-1 now hold
!
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
!              Perform a rank-2 update of A(1:k-2,1:k-2) as
!
!              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T
!                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T
!
               IF( K.GT.2 ) THEN
!
                  D12 = A( K-1, K )
                  D22 = A( K-1, K-1 ) / D12
                  D11 = A( K, K ) / D12
                  T = ONE / ( D11*D22-ONE )
                  D12 = T / D12
!
                  DO 30 J = K - 2, 1, -1
                     WKM1 = D12*( D11*A( J, K-1 )-A( J, K ) )
                     WK = D12*( D22*A( J, K )-A( J, K-1 ) )
                     DO 20 I = J, 1, -1
                        A( I, J ) = A( I, J ) - A( I, K )*WK - &
                             A( I, K-1 )*WKM1
   20                CONTINUE
                     A( J, K ) = WK
                     A( J, K-1 ) = WKM1
   30             CONTINUE
!
               END IF
!
            END IF
         END IF
!
!        Store details of the interchanges in IPIV
!
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
!
!        Decrease K and return to the start of the main loop
!
         K = K - KSTEP
         GO TO 10
!
      ELSE
!
!        Factorize A as L*D*L**T using the lower triangle of A
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2
!
         K = 1
   40    CONTINUE
!
!        If K > N, exit from loop
!
         IF( K.GT.N ) GO TO 70
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = ABS( A( K, K ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value.
!        Determine both COLMAX and IMAX.
!
         IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
!
         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN
!
!           Column K is zero or underflow, or contains a NaN:
!           set INFO and continue
!
            IF( INFO.EQ.0 ) INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            ELSE
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = K - 1 + IDAMAX( IMAX-K, A( IMAX, K ), LDA )
               ROWMAX = ABS( A( IMAX, JMAX ) )
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + IDAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( A( JMAX, IMAX ) ) )
               END IF
!
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               ELSE IF( ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
               ELSE
!
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
!
            KK = K + KSTEP - 1
            IF( KP.NE.KK ) THEN
!
!              Interchange rows and columns KK and KP in the trailing
!              submatrix A(k:n,k:n)
!
               IF( KP.LT.N ) &
                    CALL DSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
               CALL DSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), &
                    LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
            END IF
!
!           Update the trailing submatrix
!
            IF( KSTEP.EQ.1 ) THEN
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = L(k)*D(k)
!
!              where L(k) is the k-th column of L
!
               IF( K.LT.N ) THEN
!
!                 Perform a rank-1 update of A(k+1:n,k+1:n) as
!
!                 A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T
!
                  D11 = ONE / A( K, K )
                  CALL DSYR( UPLO, N-K, -D11, A( K+1, K ), 1, &
                       A( K+1, K+1 ), LDA )
!
!                 Store L(k) in column K
!
                  CALL DSCAL( N-K, D11, A( K+1, K ), 1 )
               END IF
            ELSE
!
!              2-by-2 pivot block D(k)
!
               IF( K.LT.N-1 ) THEN
!
!                 Perform a rank-2 update of A(k+2:n,k+2:n) as
!
!                 A := A - ( (A(k) A(k+1))*D(k)**(-1) ) * (A(k) A(k+1))**T
!
!                 where L(k) and L(k+1) are the k-th and (k+1)-th
!                 columns of L
!
                  D21 = A( K+1, K )
                  D11 = A( K+1, K+1 ) / D21
                  D22 = A( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
!
                  DO 60 J = K + 2, N
!
                     WK = D21*( D11*A( J, K )-A( J, K+1 ) )
                     WKP1 = D21*( D22*A( J, K+1 )-A( J, K ) )
!
                     DO 50 I = J, N
                        A( I, J ) = A( I, J ) - A( I, K )*WK - &
                             A( I, K+1 )*WKP1
   50                CONTINUE
!
                     A( J, K ) = WK
                     A( J, K+1 ) = WKP1
!
   60             CONTINUE
               END IF
            END IF
         END IF
!
!        Store details of the interchanges in IPIV
!
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
!
!        Increase K and return to the start of the main loop
!
         K = K + KSTEP
         GO TO 40
!
      END IF
!
   70 CONTINUE
!
      RETURN
!
!     End of DSYTF2
!
      END
!> \brief \b DSYTRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DSYTRF + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrf.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrf.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrf.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYTRF computes the factorization of a real symmetric matrix A using
!> the Bunch-Kaufman diagonal pivoting method.  The form of the
!> factorization is
!>
!>    A = U*D*U**T  or  A = L*D*L**T
!>
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, and D is symmetric and block diagonal with
!> 1-by-1 and 2-by-2 diagonal blocks.
!>
!> This is the blocked version of the algorithm, calling Level 3 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, the block diagonal matrix D and the multipliers used
!>          to obtain the factor U or L (see below for further details).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D.
!>          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>          interchanged and D(k,k) is a 1-by-1 diagonal block.
!>          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!>          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!>          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!>          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!>          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of WORK.  LWORK >=1.  For best performance
!>          LWORK >= N*NB, where NB is the block size returned by ILAENV.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
!>                has been completed, but the block diagonal matrix D is
!>                exactly singular, and division by zero will occur if it
!>                is used to solve a system of equations.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup doubleSYcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', then A = U*D*U**T, where
!>     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!>
!>             (   I    v    0   )   k-s
!>     U(k) =  (   0    I    0   )   s
!>             (   0    0    I   )   n-k
!>                k-s   s   n-k
!>
!>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!>  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!>
!>  If UPLO = 'L', then A = L*D*L**T, where
!>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!>
!>             (   I    0     0   )  k-1
!>     L(k) =  (   0    I     0   )  s
!>             (   0    v     I   )  n-k-s+1
!>                k-1   s  n-k-s+1
!>
!>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            IINFO, IWS, J, K, KB, LDWORK, LWKOPT, NB, NBMIN
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            ILAENV
!      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLASYF, DSYTF2, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size
!
         NB = ILAENV( 1, 'DSYTRF', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = LDWORK*NB
         IF( LWORK.LT.IWS ) THEN
            NB = MAX( LWORK / LDWORK, 1 )
            NBMIN = MAX( 2, ILAENV( 2, 'DSYTRF', UPLO, N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = 1
      END IF
      IF( NB.LT.NBMIN ) NB = N
!
      IF( UPPER ) THEN
!
!        Factorize A as U*D*U**T using the upper triangle of A
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        KB, where KB is the number of columns factorized by DLASYF;
!        KB is either NB or NB-1, or K for the last block
!
         K = N
   10    CONTINUE
!
!        If K < 1, exit from loop
!
         IF( K.LT.1 ) GO TO 40
!
         IF( K.GT.NB ) THEN
!
!           Factorize columns k-kb+1:k of A and use blocked code to
!           update columns 1:k-kb
!
            CALL DLASYF( UPLO, K, NB, KB, A, LDA, IPIV, WORK, LDWORK, &
                 IINFO )
         ELSE
!
!           Use unblocked code to factorize columns 1:k of A
!
            CALL DSYTF2( UPLO, K, A, LDA, IPIV, IINFO )
            KB = K
         END IF
!
!        Set INFO on the first occurrence of a zero pivot
!
         IF( INFO.EQ.0 .AND. IINFO.GT.0 ) &
              INFO = IINFO
!
!        Decrease K and return to the start of the main loop
!
         K = K - KB
         GO TO 10
!
      ELSE
!
!        Factorize A as L*D*L**T using the lower triangle of A
!
!        K is the main loop index, increasing from 1 to N in steps of
!        KB, where KB is the number of columns factorized by DLASYF;
!        KB is either NB or NB-1, or N-K+1 for the last block
!
         K = 1
   20    CONTINUE
!
!        If K > N, exit from loop
!
         IF( K.GT.N ) GO TO 40
!
         IF( K.LE.N-NB ) THEN
!
!           Factorize columns k:k+kb-1 of A and use blocked code to
!           update columns k+kb:n
!
            CALL DLASYF( UPLO, N-K+1, NB, KB, A( K, K ), LDA, IPIV( K ), &
                 WORK, LDWORK, IINFO )
         ELSE
!
!           Use unblocked code to factorize columns k:n of A
!
            CALL DSYTF2( UPLO, N-K+1, A( K, K ), LDA, IPIV( K ), IINFO )
            KB = N - K + 1
         END IF
!
!        Set INFO on the first occurrence of a zero pivot
!
         IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + K - 1
!
!        Adjust IPIV
!
         DO 30 J = K, K + KB - 1
            IF( IPIV( J ).GT.0 ) THEN
               IPIV( J ) = IPIV( J ) + K - 1
            ELSE
               IPIV( J ) = IPIV( J ) - K + 1
            END IF
   30    CONTINUE
!
!        Increase K and return to the start of the main loop
!
         K = K + KB
         GO TO 20
!
      END IF
!
   40 CONTINUE
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DSYTRF
!
      END
!> \brief \b DSYTRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DSYTRI + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytri.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytri.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytri.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYTRI computes the inverse of a real symmetric indefinite matrix
!> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
!> DSYTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are stored
!>          as an upper or lower triangular matrix.
!>          = 'U':  Upper triangular, form is A = U*D*U**T;
!>          = 'L':  Lower triangular, form is A = L*D*L**T.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the block diagonal matrix D and the multipliers
!>          used to obtain the factor U or L as computed by DSYTRF.
!>
!>          On exit, if INFO = 0, the (symmetric) inverse of the original
!>          matrix.  If UPLO = 'U', the upper triangular part of the
!>          inverse is formed and the part of A below the diagonal is not
!>          referenced; if UPLO = 'L' the lower triangular part of the
!>          inverse is formed and the part of A above the diagonal is
!>          not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by DSYTRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
!>               inverse could not be computed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup doubleSYcomputational
!
!  =====================================================================
      SUBROUTINE DSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K, KP, KSTEP
      DOUBLE PRECISION   AK, AKKP1, AKP1, D, T, TEMP
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      DOUBLE PRECISION   DDOT
!      EXTERNAL           LSAME, DDOT
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DSWAP, DSYMV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRI', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
!     Check that the diagonal matrix D is nonsingular.
!
      IF( UPPER ) THEN
!
!        Upper triangular storage: examine D from bottom to top
!
         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO ) RETURN
   10    CONTINUE
      ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
         DO 20 INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO ) RETURN
   20    CONTINUE
      END IF
      INFO = 0
!
      IF( UPPER ) THEN
!
!        Compute inv(A) from the factorization A = U*D*U**T.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = 1
   30    CONTINUE
!
!        If K > N, exit from loop.
!
         IF( K.GT.N ) GO TO 40
!
         IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
            A( K, K ) = ONE / A( K, K )
!
!           Compute column K of the inverse.
!
            IF( K.GT.1 ) THEN
               CALL DCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, &
                    A( 1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( K-1, WORK, 1, A( 1, K ), &
                    1 )
            END IF
            KSTEP = 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
            T = ABS( A( K, K+1 ) )
            AK = A( K, K ) / T
            AKP1 = A( K+1, K+1 ) / T
            AKKP1 = A( K, K+1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K, K ) = AKP1 / D
            A( K+1, K+1 ) = AK / D
            A( K, K+1 ) = -AKKP1 / D
!
!           Compute columns K and K+1 of the inverse.
!
            IF( K.GT.1 ) THEN
               CALL DCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, &
                    A( 1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( K-1, WORK, 1, A( 1, K ), &
                    1 )
               A( K, K+1 ) = A( K, K+1 ) -&
                    DDOT( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 )
               CALL DCOPY( K-1, A( 1, K+1 ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, &
                    A( 1, K+1 ), 1 )
               A( K+1, K+1 ) = A( K+1, K+1 ) -&
                    DDOT( K-1, WORK, 1, A( 1, K+1 ), 1 )
            END IF
            KSTEP = 2
         END IF
!
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
!
!           Interchange rows and columns K and KP in the leading
!           submatrix A(1:k+1,1:k+1)
!
            CALL DSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )
            CALL DSWAP( K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = A( K, K+1 )
               A( K, K+1 ) = A( KP, K+1 )
               A( KP, K+1 ) = TEMP
            END IF
         END IF
!
         K = K + KSTEP
         GO TO 30
   40    CONTINUE
!
      ELSE
!
!        Compute inv(A) from the factorization A = L*D*L**T.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = N
   50    CONTINUE
!
!        If K < 1, exit from loop.
!
         IF( K.LT.1 ) GO TO 60
!
         IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
            A( K, K ) = ONE / A( K, K )
!
!           Compute column K of the inverse.
!
            IF( K.LT.N ) THEN
               CALL DCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, &
                    ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( N-K, WORK, 1, A( K+1, K ), &
                    1 )
            END IF
            KSTEP = 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
            T = ABS( A( K, K-1 ) )
            AK = A( K-1, K-1 ) / T
            AKP1 = A( K, K ) / T
            AKKP1 = A( K, K-1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K-1, K-1 ) = AKP1 / D
            A( K, K ) = AK / D
            A( K, K-1 ) = -AKKP1 / D
!
!           Compute columns K-1 and K of the inverse.
!
            IF( K.LT.N ) THEN
               CALL DCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, &
                    ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( N-K, WORK, 1, A( K+1, K ), &
                    1 )
               A( K, K-1 ) = A( K, K-1 ) -&
                    DDOT( N-K, A( K+1, K ), 1, A( K+1, K-1 ), &
                    1 )
               CALL DCOPY( N-K, A( K+1, K-1 ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, &
                    ZERO, A( K+1, K-1 ), 1 )
               A( K-1, K-1 ) = A( K-1, K-1 ) - &
                    DDOT( N-K, WORK, 1, A( K+1, K-1 ), 1 )
            END IF
            KSTEP = 2
         END IF
!
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
!
!           Interchange rows and columns K and KP in the trailing
!           submatrix A(k-1:n,k-1:n)
!
            IF( KP.LT.N ) &
                 CALL DSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )
            CALL DSWAP( KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = A( K, K-1 )
               A( K, K-1 ) = A( KP, K-1 )
               A( KP, K-1 ) = TEMP
            END IF
         END IF
!
         K = K - KSTEP
         GO TO 50
   60    CONTINUE
      END IF
!
      RETURN
!
!     End of DSYTRI
!
      END
!> \brief \b DTRSM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA
!       INTEGER LDA,LDB,M,N
!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),B(LDB,*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTRSM  solves one of the matrix equations
!>
!>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!>
!> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!>
!>    op( A ) = A   or   op( A ) = A**T.
!>
!> The matrix X is overwritten on B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>           On entry, SIDE specifies whether op( A ) appears on the left
!>           or right of X as follows:
!>
!>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!>
!>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the matrix A is an upper or
!>           lower triangular matrix as follows:
!>
!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!>
!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>           On entry, TRANSA specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA = 'N' or 'n'   op( A ) = A.
!>
!>              TRANSA = 'T' or 't'   op( A ) = A**T.
!>
!>              TRANSA = 'C' or 'c'   op( A ) = A**T.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not A is unit triangular
!>           as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of B. M must be at
!>           least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of B.  N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!>           zero then  A is not referenced and  B need not be set before
!>           entry.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array of DIMENSION ( LDA, k ),
!>           where k is m when SIDE = 'L' or 'l'  
!>             and k is n when SIDE = 'R' or 'r'.
!>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!>           upper triangular part of the array  A must contain the upper
!>           triangular matrix  and the strictly lower triangular part of
!>           A is not referenced.
!>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!>           lower triangular part of the array  A must contain the lower
!>           triangular matrix  and the strictly upper triangular part of
!>           A is not referenced.
!>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!>           A  are not referenced either,  but are assumed to be  unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!>           then LDA must be at least max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!>           Before entry,  the leading  m by n part of the array  B must
!>           contain  the  right-hand  side  matrix  B,  and  on exit  is
!>           overwritten by the solution matrix  X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
!>           max( 1, m ).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level3
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 3 Blas routine.
!>
!>
!>  -- Written on 8-February-1989.
!>     Jack Dongarra, Argonne National Laboratory.
!>     Iain Duff, AERE Harwell.
!>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!>     Sven Hammarling, Numerical Algorithms Group Ltd.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!  -- Reference BLAS level3 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
!     ..
!
!  =====================================================================
!
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC MAX
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!
!     Test the input parameters.
!
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
!
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
           (.NOT.LSAME(TRANSA,'T')) .AND. &
           (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRSM ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
!
!     And when  alpha.eq.zero.
!
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
!
!     Start the operations.
!
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*inv( A )*B.
!
              IF (UPPER) THEN
                  DO 60 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      END IF
                      DO 50 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      END IF
                      DO 90 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
!
!           Form  B := alpha*inv( A**T )*B.
!
              IF (UPPER) THEN
                  DO 130 J = 1,N
                      DO 120 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          DO 110 K = 1,I - 1
                              TEMP = TEMP - A(K,I)*B(K,J)
  110                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  120                 CONTINUE
  130             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          DO 140 K = I + 1,M
                              TEMP = TEMP - A(K,I)*B(K,J)
  140                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*B*inv( A ).
!
              IF (UPPER) THEN
                  DO 210 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 170 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  170                     CONTINUE
                      END IF
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 200 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  200                     CONTINUE
                      END IF
  210             CONTINUE
              ELSE
                  DO 260 J = N,1,-1
                      IF (ALPHA.NE.ONE) THEN
                          DO 220 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  220                     CONTINUE
                      END IF
                      DO 240 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              DO 230 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 250 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              END IF
          ELSE
!
!           Form  B := alpha*B*inv( A**T ).
!
              IF (UPPER) THEN
                  DO 310 K = N,1,-1
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
                      DO 290 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 280 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  280                         CONTINUE
                          END IF
  290                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 300 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  300                     CONTINUE
                      END IF
  310             CONTINUE
              ELSE
                  DO 360 K = 1,N
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 320 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  320                     CONTINUE
                      END IF
                      DO 340 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 330 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  330                         CONTINUE
                          END IF
  340                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 350 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  350                     CONTINUE
                      END IF
  360             CONTINUE
              END IF
          END IF
      END IF
!
      RETURN
!
!     End of DTRSM .
!
      END
!> \brief \b IDAMAX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION IDAMAX(N,DX,INCX)
! 
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    IDAMAX finds the index of the first element having maximum absolute value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2015
!
!> \ingroup aux_blas
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
!
!  -- Reference BLAS level1 routine (version 3.6.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DABS
!     ..
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
!
!        code for increment equal to 1
!
         DMAX = DABS(DX(1))
         DO I = 2,N
            IF (DABS(DX(I)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(I))
            END IF
         END DO
      ELSE
!
!        code for increment not equal to 1
!
         IX = 1
         DMAX = DABS(DX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (DABS(DX(IX)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END
!> \brief \b IEEECK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download IEEECK + dependencies 
!> [TGZ]</a> 
!> [ZIP]</a> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
! 
!       .. Scalar Arguments ..
!       INTEGER            ISPEC
!       REAL               ONE, ZERO
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> IEEECK is called from the ILAENV to verify that Infinity and
!> possibly NaN arithmetic is safe (i.e. will not trap).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          Specifies whether to test just for inifinity arithmetic
!>          or whether to test for infinity and NaN arithmetic.
!>          = 0: Verify infinity arithmetic only.
!>          = 1: Verify infinity and NaN arithmetic.
!> \endverbatim
!>
!> \param[in] ZERO
!> \verbatim
!>          ZERO is REAL
!>          Must contain the value 0.0
!>          This is passed to prevent the compiler from optimizing
!>          away this code.
!> \endverbatim
!>
!> \param[in] ONE
!> \verbatim
!>          ONE is REAL
!>          Must contain the value 1.0
!>          This is passed to prevent the compiler from optimizing
!>          away this code.
!>
!>  RETURN VALUE:  INTEGER
!>          = 0:  Arithmetic failed to produce the correct answers
!>          = 1:  Arithmetic produced the correct answers
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, &
           NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
      IEEECK = 1
!
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      IF( ISPEC.EQ.0 ) RETURN
!
      NAN1 = POSINF + NEGINF
!
      NAN2 = POSINF / NEGINF
!
      NAN3 = POSINF / POSINF
!
      NAN4 = POSINF*ZERO
!
      NAN5 = NEGINF*NEGZRO
!
      NAN6 = NAN5*ZERO
!
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      RETURN
      END
!> \brief \b ILAENV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ILAENV + dependencies 
!> [TGZ]</a> 
!> [ZIP]</a> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
! 
!       .. Scalar Arguments ..
!       CHARACTER*( * )    NAME, OPTS
!       INTEGER            ISPEC, N1, N2, N3, N4
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILAENV is called from the LAPACK routines to choose problem-dependent
!> parameters for the local environment.  See ISPEC for a description of
!> the parameters.
!>
!> ILAENV returns an INTEGER
!> if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
!> if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!>
!> This version provides a set of parameters which should give good,
!> but not optimal, performance on many of the currently available
!> computers.  Users are encouraged to modify this subroutine to set
!> the tuning parameters for their particular machine using the option
!> and problem size information in the arguments.
!>
!> This routine will not function correctly if it is converted to all
!> lower case.  Converting it to all upper case is allowed.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          Specifies the parameter to be returned as the value of
!>          ILAENV.
!>          = 1: the optimal blocksize; if this value is 1, an unblocked
!>               algorithm will give the best performance.
!>          = 2: the minimum block size for which the block routine
!>               should be used; if the usable block size is less than
!>               this value, an unblocked routine should be used.
!>          = 3: the crossover point (in a block routine, for N less
!>               than this value, an unblocked routine should be used)
!>          = 4: the number of shifts, used in the nonsymmetric
!>               eigenvalue routines (DEPRECATED)
!>          = 5: the minimum column dimension for blocking to be used;
!>               rectangular blocks must have dimension at least k by m,
!>               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!>          = 6: the crossover point for the SVD (when reducing an m by n
!>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!>               this value, a QR factorization is used first to reduce
!>               the matrix to a triangular form.)
!>          = 7: the number of processors
!>          = 8: the crossover point for the multishift QR method
!>               for nonsymmetric eigenvalue problems (DEPRECATED)
!>          = 9: maximum size of the subproblems at the bottom of the
!>               computation tree in the divide-and-conquer algorithm
!>               (used by xGELSD and xGESDD)
!>          =10: ieee NaN arithmetic can be trusted not to trap
!>          =11: infinity arithmetic can be trusted not to trap
!>          12 <= ISPEC <= 16:
!>               xHSEQR or related subroutines,
!>               see IPARMQ for detailed explanation
!> \endverbatim
!>
!> \param[in] NAME
!> \verbatim
!>          NAME is CHARACTER*(*)
!>          The name of the calling subroutine, in either upper case or
!>          lower case.
!> \endverbatim
!>
!> \param[in] OPTS
!> \verbatim
!>          OPTS is CHARACTER*(*)
!>          The character options to the subroutine NAME, concatenated
!>          into a single character string.  For example, UPLO = 'U',
!>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!>          be specified as OPTS = 'UTN'.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!> \endverbatim
!>
!> \param[in] N3
!> \verbatim
!>          N3 is INTEGER
!> \endverbatim
!>
!> \param[in] N4
!> \verbatim
!>          N4 is INTEGER
!>          Problem dimensions for the subroutine NAME; these may not all
!>          be required.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2015
!
!> \ingroup auxOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The following conventions have been used when calling ILAENV from the
!>  LAPACK routines:
!>  1)  OPTS is a concatenation of all of the character options to
!>      subroutine NAME, in the same order that they appear in the
!>      argument list for NAME, even if they are not used in determining
!>      the value of the parameter specified by ISPEC.
!>  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!>      that they appear in the argument list for NAME.  N1 is used
!>      first, N2 second, and so on, and unused problem dimensions are
!>      passed a value of -1.
!>  3)  The parameter value returned by ILAENV is checked for validity in
!>      the calling subroutine.  For example, ILAENV is used to retrieve
!>      the optimal blocksize for STRTRI as follows:
!>
!>      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!>      IF( NB.LE.1 ) NB = MAX( 1, N )
!> \endverbatim
!>
!  =====================================================================
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!  -- LAPACK auxiliary routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. External Functions ..
!      INTEGER            IEEECK, IPARMQ
!      EXTERNAL           IEEECK, IPARMQ
!     ..
!     .. Executable Statements ..
!
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120, &
           130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
   10 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) &
                    SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
              ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
              ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                    ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                    ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I: &
                    I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 ) &
                    SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) ) RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
!
      GO TO ( 50, 60, 70 )ISPEC
!
   50 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
              C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      ELSE IF( C2.EQ.'GG' ) THEN
         NB = 32
         IF( C3.EQ.'HD3' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         END IF
      END IF
      ILAENV = NB
      RETURN
!
   60 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ. &
              'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'GG' ) THEN
         NBMIN = 2
         IF( C3.EQ.'HD3' ) THEN
            NBMIN = 2
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
   70 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ. &
              'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'GG' ) THEN
         NX = 128
         IF( C3.EQ.'HD3' ) THEN
            NX = 128
         END IF
      END IF
      ILAENV = NX
      RETURN
!
   80 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
   90 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  100 CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  110 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  120 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
  130 CONTINUE
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
      ILAENV = 25
      RETURN
!
  140 CONTINUE
!
!     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
!
  150 CONTINUE
!
!     ISPEC = 11: infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
!
  160 CONTINUE
!
!     12 <= ISPEC <= 16: xHSEQR or related subroutines.
!
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
!
!     End of ILAENV
!
   END FUNCTION ILAENV
!> \brief \b IPARMQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download IPARMQ + dependencies 
!> [TGZ]</a> 
!> [ZIP]</a> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
! 
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, ISPEC, LWORK, N
!       CHARACTER          NAME*( * ), OPTS*( * )
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      This program sets problem and machine dependent parameters
!>      useful for xHSEQR and related subroutines for eigenvalue
!>      problems. It is called whenever
!>      IPARMQ is called with 12 <= ISPEC <= 16
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is integer scalar
!>              ISPEC specifies which tunable parameter IPARMQ should
!>              return.
!>
!>              ISPEC=12: (INMIN)  Matrices of order nmin or less
!>                        are sent directly to xLAHQR, the implicit
!>                        double shift QR algorithm.  NMIN must be
!>                        at least 11.
!>
!>              ISPEC=13: (INWIN)  Size of the deflation window.
!>                        This is best set greater than or equal to
!>                        the number of simultaneous shifts NS.
!>                        Larger matrices benefit from larger deflation
!>                        windows.
!>
!>              ISPEC=14: (INIBL) Determines when to stop nibbling and
!>                        invest in an (expensive) multi-shift QR sweep.
!>                        If the aggressive early deflation subroutine
!>                        finds LD converged eigenvalues from an order
!>                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
!>                        then the next QR sweep is skipped and early
!>                        deflation is applied immediately to the
!>                        remaining active diagonal block.  Setting
!>                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
!>                        multi-shift QR sweep whenever early deflation
!>                        finds a converged eigenvalue.  Setting
!>                        IPARMQ(ISPEC=14) greater than or equal to 100
!>                        prevents TTQRE from skipping a multi-shift
!>                        QR sweep.
!>
!>              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
!>                        a multi-shift QR iteration.
!>
!>              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
!>                        following meanings.
!>                        0:  During the multi-shift QR/QZ sweep,
!>                            blocked eigenvalue reordering, blocked
!>                            Hessenberg-triangular reduction,
!>                            reflections and/or rotations are not
!>                            accumulated when updating the
!>                            far-from-diagonal matrix entries.
!>                        1:  During the multi-shift QR/QZ sweep,
!>                            blocked eigenvalue reordering, blocked
!>                            Hessenberg-triangular reduction,
!>                            reflections and/or rotations are
!>                            accumulated, and matrix-matrix
!>                            multiplication is used to update the
!>                            far-from-diagonal matrix entries.
!>                        2:  During the multi-shift QR/QZ sweep,
!>                            blocked eigenvalue reordering, blocked
!>                            Hessenberg-triangular reduction,
!>                            reflections and/or rotations are
!>                            accumulated, and 2-by-2 block structure
!>                            is exploited during matrix-matrix
!>                            multiplies.
!>                        (If xTRMM is slower than xGEMM, then
!>                        IPARMQ(ISPEC=16)=1 may be more efficient than
!>                        IPARMQ(ISPEC=16)=2 despite the greater level of
!>                        arithmetic work implied by the latter choice.)
!> \endverbatim
!>
!> \param[in] NAME
!> \verbatim
!>          NAME is character string
!>               Name of the calling subroutine
!> \endverbatim
!>
!> \param[in] OPTS
!> \verbatim
!>          OPTS is character string
!>               This is a concatenation of the string arguments to
!>               TTQRE.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer scalar
!>               N is the order of the Hessenberg matrix H.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>               It is assumed that H is already upper triangular
!>               in rows and columns 1:ILO-1 and IHI+1:N.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is integer scalar
!>               The amount of workspace available.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2015
!
!> \ingroup auxOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>       Little is known about how best to choose these parameters.
!>       It is possible to use different values of the parameters
!>       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
!>
!>       It is probably best to choose different parameters for
!>       different matrices and different parameters at different
!>       times during the iteration, but this has not been
!>       implemented --- yet.
!>
!>
!>       The best choices of most of the parameters depend
!>       in an ill-understood way on the relative execution
!>       rate of xLAQR3 and xLAQR5 and on the nature of each
!>       particular eigenvalue problem.  Experiment may be the
!>       only practical way to determine which choices are most
!>       effective.
!>
!>       Following is a list of default values supplied by IPARMQ.
!>       These defaults may be adjusted in order to attain better
!>       performance in any particular computational environment.
!>
!>       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
!>                        Default: 75. (Must be at least 11.)
!>
!>       IPARMQ(ISPEC=13) Recommended deflation window size.
!>                        This depends on ILO, IHI and NS, the
!>                        number of simultaneous shifts returned
!>                        by IPARMQ(ISPEC=15).  The default for
!>                        (IHI-ILO+1).LE.500 is NS.  The default
!>                        for (IHI-ILO+1).GT.500 is 3*NS/2.
!>
!>       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
!>
!>       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
!>                        a multi-shift QR iteration.
!>
!>                        If IHI-ILO+1 is ...
!>
!>                        greater than      ...but less    ... the
!>                        or equal to ...      than        default is
!>
!>                                0               30       NS =   2+
!>                               30               60       NS =   4+
!>                               60              150       NS =  10
!>                              150              590       NS =  **
!>                              590             3000       NS =  64
!>                             3000             6000       NS = 128
!>                             6000             infinity   NS = 256
!>
!>                    (+)  By default matrices of this order are
!>                         passed to the implicit double shift routine
!>                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
!>                         values of NS are used only in case of a rare
!>                         xLAHQR failure.
!>
!>                    (**) The asterisks (**) indicate an ad-hoc
!>                         function increasing from 10 to 64.
!>
!>       IPARMQ(ISPEC=16) Select structured matrix multiply.
!>                        (See ISPEC=16 above for details.)
!>                        Default: 3.
!> \endverbatim
!>
!  =====================================================================
      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
!
!  -- LAPACK auxiliary routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
!
!  ================================================================
!     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14, &
           ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14, &
           NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
!     ..
!     .. Local Scalars ..
      INTEGER            NH, NS
      INTEGER            I, IC, IZ
      CHARACTER          SUBNAM*6
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
!     ..
!     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR. &
           ( ISPEC.EQ.IACC22 ) ) THEN
!
!        ==== Set the number simultaneous shifts ====
!
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 ) NS = 4
         IF( NH.GE.60 ) NS = 10
         IF( NH.GE.150 ) &
              NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
         IF( NH.GE.590 ) NS = 64
         IF( NH.GE.3000 ) NS = 128
         IF( NH.GE.6000 ) NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
!
      IF( ISPEC.EQ.INMIN ) THEN
!
!
!        ===== Matrices of order smaller than NMIN get sent
!        .     to xLAHQR, the classic double shift algorithm.
!        .     This must be at least 11. ====
!
         IPARMQ = NMIN
!
      ELSE IF( ISPEC.EQ.INIBL ) THEN
!
!        ==== INIBL: skip a multi-shift qr iteration and
!        .    whenever aggressive early deflation finds
!        .    at least (NIBBLE*(window size)/100) deflations. ====
!
         IPARMQ = NIBBLE
!
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
!
!        ==== NSHFTS: The number of simultaneous shifts =====
!
         IPARMQ = NS
!
      ELSE IF( ISPEC.EQ.INWIN ) THEN
!
!        ==== NW: deflation window size.  ====
!
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
!
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
!
!        ==== IACC22: Whether to accumulate reflections
!        .     before updating the far-from-diagonal elements
!        .     and whether to use 2-by-2 block structure while
!        .     doing it.  A small amount of work could be saved
!        .     by making this choice dependent also upon the
!        .     NH=IHI-ILO+1.
!
!
!        Convert NAME to upper case if the first character is lower case.
!
         IPARMQ = 0
         SUBNAM = NAME
         IC = ICHAR( SUBNAM( 1: 1 ) )
         IZ = ICHAR( 'Z' )
         IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!           ASCII character set
!
            IF( IC.GE.97 .AND. IC.LE.122 ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( IC.GE.97 .AND. IC.LE.122 ) &
                       SUBNAM( I: I ) = CHAR( IC-32 )
               END DO
            END IF
!
         ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!           EBCDIC character set
!
            IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                 ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                 ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC+64 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                       ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                       ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I: &
                       I ) = CHAR( IC+64 )
               END DO
            END IF
!
         ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!           Prime machines:  ASCII+128
!
            IF( IC.GE.225 .AND. IC.LE.250 ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( IC.GE.225 .AND. IC.LE.250 ) &
                       SUBNAM( I: I ) = CHAR( IC-32 )
               END DO
            END IF
         END IF
!
         IF( SUBNAM( 2:6 ).EQ.'GGHRD' .OR. &
              SUBNAM( 2:6 ).EQ.'GGHD3' ) THEN
            IPARMQ = 1
            IF( NH.GE.K22MIN ) IPARMQ = 2
         ELSE IF ( SUBNAM( 4:6 ).EQ.'EXC' ) THEN
            IF( NH.GE.KACMIN ) IPARMQ = 1
            IF( NH.GE.K22MIN ) IPARMQ = 2
         ELSE IF ( SUBNAM( 2:6 ).EQ.'HSEQR' .OR. &
              SUBNAM( 2:5 ).EQ.'LAQR' ) THEN
            IF( NS.GE.KACMIN ) IPARMQ = 1
            IF( NS.GE.K22MIN ) IPARMQ = 2
         END IF
!
      ELSE
!        ===== invalid value of ispec =====
         IPARMQ = -1
!
      END IF
!
!     ==== End of IPARMQ ====
!
      END
!> \brief \b LSAME
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION LSAME(CA,CB)
! 
!       .. Scalar Arguments ..
!       CHARACTER CA,CB
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> LSAME returns .TRUE. if CA is the same letter as CB regardless of
!> case.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CA
!> \verbatim
!>          CA is CHARACTER*1
!> \endverbatim
!>
!> \param[in] CB
!> \verbatim
!>          CB is CHARACTER*1
!>          CA and CB specify the single characters to be compared.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup aux_blas
!
!  =====================================================================
      LOGICAL FUNCTION LSAME(CA,CB)
!
!  -- Reference BLAS level1 routine (version 3.1) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER CA,CB
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER INTA,INTB,ZCODE
!     ..
!
!     Test if the characters are equal
!
      LSAME = CA .EQ. CB
      IF (LSAME) RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR('Z')
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
!
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
          IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
          IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
!
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
          IF (INTA.GE.129 .AND. INTA.LE.137 .OR. &
               INTA.GE.145 .AND. INTA.LE.153 .OR. &
               INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
          IF (INTB.GE.129 .AND. INTB.LE.137 .OR. &
               INTB.GE.145 .AND. INTB.LE.153 .OR. &
               INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
!
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
          IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
          IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
      LSAME = INTA .EQ. INTB
!
!     RETURN
!
!     End of LSAME
!
      END
!> \brief \b LSAMEN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download LSAMEN + dependencies 
!> [TGZ]</a> 
!> [ZIP]</a> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION LSAMEN( N, CA, CB )
! 
!       .. Scalar Arguments ..
!       CHARACTER*( * )    CA, CB
!       INTEGER            N
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> LSAMEN  tests if the first N letters of CA are the same as the
!> first N letters of CB, regardless of case.
!> LSAMEN returns .TRUE. if CA and CB are equivalent except for case
!> and .FALSE. otherwise.  LSAMEN also returns .FALSE. if LEN( CA )
!> or LEN( CB ) is less than N.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of characters in CA and CB to be compared.
!> \endverbatim
!>
!> \param[in] CA
!> \verbatim
!>          CA is CHARACTER*(*)
!> \endverbatim
!>
!> \param[in] CB
!> \verbatim
!>          CB is CHARACTER*(*)
!>          CA and CB specify two character strings of length at least N.
!>          Only the first N characters of each string will be accessed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      LOGICAL          FUNCTION LSAMEN( N, CA, CB )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    CA, CB
      INTEGER            N
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          LEN
!     ..
!     .. Executable Statements ..
!
      LSAMEN = .FALSE.
      IF( LEN( CA ).LT.N .OR. LEN( CB ).LT.N ) GO TO 20
!
!     Do for each character in the two strings.
!
      DO 10 I = 1, N
!
!        Test if the characters are equal using LSAME.
!
         IF( .NOT.LSAME( CA( I: I ), CB( I: I ) ) ) GO TO 20
!
   10 CONTINUE
      LSAMEN = .TRUE.
!
   20 CONTINUE
      RETURN
!
!     End of LSAMEN
!
      END
!> \brief \b XERBLA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE XERBLA( SRNAME, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER*(*)      SRNAME
!       INTEGER            INFO
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> XERBLA  is an error handler for the LAPACK routines.
!> It is called by an LAPACK routine if an input parameter has an
!> invalid value.  A message is printed and execution stops.
!>
!> Installers may consider modifying the STOP statement in order to
!> call system-specific exception-handling facilities.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SRNAME
!> \verbatim
!>          SRNAME is CHARACTER*(*)
!>          The name of the routine which called XERBLA.
!> \endverbatim
!>
!> \param[in] INFO
!> \verbatim
!>          INFO is INTEGER
!>          The position of the invalid parameter in the parameter list
!>          of the calling routine.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup aux_blas
!
!  =====================================================================
      SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          LEN_TRIM
!     ..
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ', &
           'an illegal value' )
!
!     End of XERBLA
!
    END SUBROUTINE XERBLA
!
!> \brief <b> DSPEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DSPEVD + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dspevd.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dspevd.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dspevd.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK,
!                          IWORK, LIWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, LDZ, LIWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSPEVD computes all the eigenvalues and, optionally, eigenvectors
!> of a real symmetric matrix A in packed storage. If eigenvectors are
!> desired, it uses a divide and conquer algorithm.
!>
!> The divide and conquer algorithm makes very mild assumptions about
!> floating point arithmetic. It will work on machines with a guard
!> digit in add/subtract, or on those binary machines without guard
!> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!> Cray-2. It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBZ
!> \verbatim
!>          JOBZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only;
!>          = 'V':  Compute eigenvalues and eigenvectors.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the symmetric matrix
!>          A, packed columnwise in a linear array.  The j-th column of A
!>          is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!>
!>          On exit, AP is overwritten by values generated during the
!>          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
!>          and first superdiagonal of the tridiagonal matrix T overwrite
!>          the corresponding elements of A, and if UPLO = 'L', the
!>          diagonal and first subdiagonal of T overwrite the
!>          corresponding elements of A.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
!>          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
!>          eigenvectors of the matrix A, with the i-th column of Z
!>          holding the eigenvector associated with W(i).
!>          If JOBZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          JOBZ = 'V', LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array,
!>                                         dimension (LWORK)
!>          On exit, if INFO = 0, WORK(1) returns the required LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If N <= 1,               LWORK must be at least 1.
!>          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N.
!>          If JOBZ = 'V' and N > 1, LWORK must be at least
!>                                                 1 + 6*N + N**2.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the required sizes of the WORK and IWORK
!>          arrays, returns these values as the first entries of the WORK
!>          and IWORK arrays, and no error message related to LWORK or
!>          LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the required LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.
!>          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the required sizes of the WORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK and IWORK arrays, and no error message related to
!>          LWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = i, the algorithm failed to converge; i
!>                off-diagonal elements of an intermediate tridiagonal
!>                form did not converge to zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup doubleOTHEReigen
!
!  =====================================================================
  SUBROUTINE DSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, &
       IWORK, LIWORK, INFO )
!
!  -- LAPACK driver routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, LDZ, LIWORK, LWORK, N
!     ..
!     .. Array Arguments ..
    INTEGER            IWORK( * )
    DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
    DOUBLE PRECISION   ZERO, ONE
    PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
    LOGICAL            LQUERY, WANTZ
    INTEGER            IINFO, INDE, INDTAU, INDWRK, ISCALE, LIWMIN,&
         LLWORK, LWMIN
    DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,&
         SMLNUM
!     ..
!     .. External Functions ..
!?    LOGICAL            LSAME
!?    DOUBLE PRECISION   DLAMCH, DLANSP
!?    EXTERNAL           LSAME, DLAMCH, DLANSP
!     ..
!     .. External Subroutines ..
!?    EXTERNAL           DOPMTR, DSCAL, DSPTRD, DSTEDC, DSTERF, XERBLA
!    EXTERNAL           DSTERF
!     ..
!     .. Intrinsic Functions ..
!?      INTRINSIC          SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
    WANTZ = LSAME( JOBZ, 'V' )
    LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
!
    INFO = 0
    IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
       INFO = -1
    ELSE IF( .NOT.( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) )&
         THEN
       INFO = -2
    ELSE IF( N.LT.0 ) THEN
       INFO = -3
    ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
       INFO = -7
    END IF
!
    IF( INFO.EQ.0 ) THEN
       IF( N.LE.1 ) THEN
          LIWMIN = 1
          LWMIN = 1
       ELSE
          IF( WANTZ ) THEN
             LIWMIN = 3 + 5*N
             LWMIN = 1 + 6*N + N**2
          ELSE
             LIWMIN = 1
             LWMIN = 2*N
          END IF
       END IF
       IWORK( 1 ) = LIWMIN
       WORK( 1 ) = LWMIN
!
       IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
          INFO = -9
       ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
          INFO = -11
       END IF
    END IF
!
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DSPEVD', -INFO )
       RETURN
    ELSE IF( LQUERY ) THEN
       RETURN
    END IF
!
!     Quick return if possible
!
    IF( N.EQ.0 )&
         RETURN
!
    IF( N.EQ.1 ) THEN
       W( 1 ) = AP( 1 )
       IF( WANTZ )&
            Z( 1, 1 ) = ONE
       RETURN
    END IF
!
!     Get machine constants.
!
    SAFMIN = DLAMCH( 'Safe minimum' )
    EPS = DLAMCH( 'Precision' )
    SMLNUM = SAFMIN / EPS
    BIGNUM = ONE / SMLNUM
    RMIN = SQRT( SMLNUM )
    RMAX = SQRT( BIGNUM )
!
!     Scale matrix to allowable range, if necessary.
!
    ANRM = DLANSP( 'M', UPLO, N, AP, WORK )
    ISCALE = 0
    IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
       ISCALE = 1
       SIGMA = RMIN / ANRM
    ELSE IF( ANRM.GT.RMAX ) THEN
       ISCALE = 1
       SIGMA = RMAX / ANRM
    END IF
    IF( ISCALE.EQ.1 ) THEN
       CALL DSCAL( ( N*( N+1 ) ) / 2, SIGMA, AP, 1 )
    END IF
!
!     Call DSPTRD to reduce symmetric packed matrix to tridiagonal form.
!
    INDE = 1
    INDTAU = INDE + N
    CALL DSPTRD( UPLO, N, AP, W, WORK( INDE ), WORK( INDTAU ), IINFO )
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, first call
!     DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
!     tridiagonal matrix, then call DOPMTR to multiply it by the
!     Householder transformations represented in AP.
!
    IF( .NOT.WANTZ ) THEN
       CALL DSTERF( N, W, WORK( INDE ), INFO )
    ELSE
       INDWRK = INDTAU + N
       LLWORK = LWORK - INDWRK + 1
       CALL DSTEDC( 'I', N, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ),&
            LLWORK, IWORK, LIWORK, INFO )
       CALL DOPMTR( 'L', UPLO, 'N', N, N, AP, WORK( INDTAU ), Z, LDZ,&
            WORK( INDWRK ), IINFO )
    END IF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
    IF( ISCALE.EQ.1 ) &
         CALL DSCAL( N, ONE / SIGMA, W, 1 )
!
    WORK( 1 ) = LWMIN
    IWORK( 1 ) = LIWMIN
    RETURN
!
!     End of DSPEVD
!
  END SUBROUTINE DSPEVD
!
!
!> \brief \b DLANSP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a symmetric matrix supplied in packed form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLANSP + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlansp.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlansp.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlansp.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLANSP( NORM, UPLO, N, AP, WORK )
! 
!       .. Scalar Arguments ..
!       CHARACTER          NORM, UPLO
!       INTEGER            N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AP( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLANSP  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> real symmetric matrix A,  supplied in packed form.
!> \endverbatim
!>
!> \return DLANSP
!> \verbatim
!>
!>    DLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!>             (
!>             ( norm1(A),         NORM = '1', 'O' or 'o'
!>             (
!>             ( normI(A),         NORM = 'I' or 'i'
!>             (
!>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normF  denotes the  Frobenius norm of a matrix (square root of sum of
!> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies the value to be returned in DLANSP as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is supplied.
!>          = 'U':  Upper triangular part of A is supplied
!>          = 'L':  Lower triangular part of A is supplied
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, DLANSP is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          The upper or lower triangle of the symmetric matrix A, packed
!>          columnwise in a linear array.  The j-th column of A is stored
!>          in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!>          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!>          WORK is not referenced.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup doubleOTHERauxiliary
!
!  =====================================================================
  DOUBLE PRECISION FUNCTION DLANSP( NORM, UPLO, N, AP, WORK )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
    CHARACTER          NORM, UPLO
    INTEGER            N
!     ..
!     .. Array Arguments ..
    DOUBLE PRECISION   AP( * ), WORK( * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
    DOUBLE PRECISION   ONE, ZERO
    PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
    INTEGER            I, J, K
    DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
!    EXTERNAL           DLASSQ
!     ..
!     .. External Functions ..
!    LOGICAL            LSAME, DISNAN
!    EXTERNAL           LSAME, DISNAN
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
    IF( N.EQ.0 ) THEN
       VALUE = ZERO
    ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
       VALUE = ZERO
       IF( LSAME( UPLO, 'U' ) ) THEN
          K = 1
          DO 20 J = 1, N
             DO 10 I = K, K + J - 1
                SUM = ABS( AP( I ) )
                IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
10           CONTINUE
             K = K + J
20        CONTINUE
       ELSE
          K = 1
          DO 40 J = 1, N
             DO 30 I = K, K + N - J
                SUM = ABS( AP( I ) )
                IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
30           CONTINUE
             K = K + N - J + 1
40        CONTINUE
       END IF
    ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ).OR.&
         ( NORM.EQ.'1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
       VALUE = ZERO
       K = 1
       IF( LSAME( UPLO, 'U' ) ) THEN
          DO 60 J = 1, N
             SUM = ZERO
             DO 50 I = 1, J - 1
                ABSA = ABS( AP( K ) )
                SUM = SUM + ABSA
                WORK( I ) = WORK( I ) + ABSA
                K = K + 1
50           CONTINUE
             WORK( J ) = SUM + ABS( AP( K ) )
             K = K + 1
60        CONTINUE
          DO 70 I = 1, N
             SUM = WORK( I )
             IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) &
                  VALUE = SUM
70        CONTINUE
       ELSE
          DO 80 I = 1, N
             WORK( I ) = ZERO
80        CONTINUE
          DO 100 J = 1, N
             SUM = WORK( J ) + ABS( AP( K ) )
             K = K + 1
             DO 90 I = J + 1, N
                ABSA = ABS( AP( K ) )
                SUM = SUM + ABSA
                WORK( I ) = WORK( I ) + ABSA
                K = K + 1
90           CONTINUE
             IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
100       CONTINUE
       END IF
    ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
       SCALE = ZERO
       SUM = ONE
       K = 2
       IF( LSAME( UPLO, 'U' ) ) THEN
          DO 110 J = 2, N
             CALL DLASSQ( J-1, AP( K ), 1,&
                  SCALE, SUM )
             K = K + J
110       CONTINUE
       ELSE
          DO 120 J = 1, N - 1
             CALL DLASSQ( N-J, AP( K ), 1, SCALE, SUM )
             K = K + N - J + 1
120       CONTINUE
       END IF
       SUM = 2*SUM
       K = 1
       DO 130 I = 1, N
          IF( AP( K ).NE.ZERO ) THEN
             ABSA = ABS( AP( K ) )
             IF( SCALE.LT.ABSA ) THEN
                SUM = ONE + SUM*( SCALE / ABSA )**2
                SCALE = ABSA
             ELSE
                SUM = SUM + ( ABSA / SCALE )**2
             END IF
          END IF
          IF( LSAME( UPLO, 'U' ) ) THEN
             K = K + I + 1
          ELSE
             K = K + N - I + 1
          END IF
130    CONTINUE
       VALUE = SCALE*SQRT( SUM )
    END IF
!
    DLANSP = VALUE
    RETURN
!
!     End of DLANSP
!
  END FUNCTION DLANSP
!
!=
!
!> \brief \b DLASSQ updates a sum of squares represented in scaled form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLASSQ + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlassq.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlassq.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlassq.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
! 
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       DOUBLE PRECISION   SCALE, SUMSQ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   X( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASSQ  returns the values  scl  and  smsq  such that
!>
!>    ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!>
!> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!> assumed to be non-negative and  scl  returns the value
!>
!>    scl = max( scale, abs( x( i ) ) ).
!>
!> scale and sumsq must be supplied in SCALE and SUMSQ and
!> scl and smsq are overwritten on SCALE and SUMSQ respectively.
!>
!> The routine makes only one pass through the vector x.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of elements to be used from the vector X.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (N)
!>          The vector for which a scaled sum of squares is computed.
!>             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of the vector X.
!>          INCX > 0.
!> \endverbatim
!>
!> \param[in,out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          On entry, the value  scale  in the equation above.
!>          On exit, SCALE is overwritten with  scl , the scaling factor
!>          for the sum of squares.
!> \endverbatim
!>
!> \param[in,out] SUMSQ
!> \verbatim
!>          SUMSQ is DOUBLE PRECISION
!>          On entry, the value  sumsq  in the equation above.
!>          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!>          squares from which  scl  has been factored out.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
  SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
    INTEGER            INCX, N
    DOUBLE PRECISION   SCALE, SUMSQ
!     ..
!     .. Array Arguments ..
    DOUBLE PRECISION   X( * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
    DOUBLE PRECISION   ZERO
    PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
    INTEGER            IX
    DOUBLE PRECISION   ABSXI
!     ..
!     .. External Functions ..
!      LOGICAL            DISNAN
!      EXTERNAL           DISNAN
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
    IF( N.GT.0 ) THEN
       DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
          ABSXI = ABS( X( IX ) )
          IF( ABSXI.GT.ZERO.OR.DISNAN( ABSXI ) ) THEN
             IF( SCALE.LT.ABSXI ) THEN
                SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                SCALE = ABSXI
             ELSE
                SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
             END IF
          END IF
10     CONTINUE
    END IF
    RETURN
!
!     End of DLASSQ
!
 END SUBROUTINE DLASSQ
!
!=
!
!> \brief \b DSPTRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DSPTRD + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsptrd.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsptrd.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsptrd.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSPTRD( UPLO, N, AP, D, E, TAU, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AP( * ), D( * ), E( * ), TAU( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSPTRD reduces a real symmetric matrix A stored in packed form to
!> symmetric tridiagonal form T by an orthogonal similarity
!> transformation: Q**T * A * Q = T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the symmetric matrix
!>          A, packed columnwise in a linear array.  The j-th column of A
!>          is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!>          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!>          of A are overwritten by the corresponding elements of the
!>          tridiagonal matrix T, and the elements above the first
!>          superdiagonal, with the array TAU, represent the orthogonal
!>          matrix Q as a product of elementary reflectors; if UPLO
!>          = 'L', the diagonal and first subdiagonal of A are over-
!>          written by the corresponding elements of the tridiagonal
!>          matrix T, and the elements below the first subdiagonal, with
!>          the array TAU, represent the orthogonal matrix Q as a product
!>          of elementary reflectors. See Further Details.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal elements of the tridiagonal matrix T:
!>          D(i) = A(i,i).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The off-diagonal elements of the tridiagonal matrix T:
!>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup doubleOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(n-1) . . . H(2) H(1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in AP,
!>  overwriting A(1:i-1,i+1), and tau is stored in TAU(i).
!>
!>  If UPLO = 'L', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(1) H(2) . . . H(n-1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,
!>  overwriting A(i+2:n,i), and tau is stored in TAU(i).
!> \endverbatim
!>
!  =====================================================================
 SUBROUTINE DSPTRD( UPLO, N, AP, D, E, TAU, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
   CHARACTER          UPLO
   INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   AP( * ), D( * ), E( * ), TAU( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE, ZERO, HALF
   PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0,&
        HALF = 1.0D0 / 2.0D0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            UPPER
   INTEGER            I, I1, I1I1, II
   DOUBLE PRECISION   ALPHA, TAUI
!     ..
!     .. External Subroutines ..
!?   EXTERNAL           DAXPY, DLARFG, DSPMV, DSPR2, XERBLA
!   EXTERNAL           DAXPY, DSPMV, DSPR2
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      DOUBLE PRECISION   DDOT
!      EXTERNAL           LSAME, DDOT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
   INFO = 0
   UPPER = LSAME( UPLO, 'U' )
   IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
   ELSE IF( N.LT.0 ) THEN
      INFO = -2
   END IF
   IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'DSPTRD', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N.LE.0 )&
        RETURN
!
   IF( UPPER ) THEN
!
!        Reduce the upper triangle of A.
!        I1 is the index in AP of A(1,I+1).
!
      I1 = N*( N-1 ) / 2 + 1
      DO 10 I = N - 1, 1, -1
!
!           Generate elementary reflector H(i) = I - tau * v * v**T
!           to annihilate A(1:i-1,i+1)
!
         CALL DLARFG( I, AP( I1+I-1 ), AP( I1 ), 1, TAUI )
         E( I ) = AP( I1+I-1 )
!
         IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
            AP( I1+I-1 ) = ONE
!
!              Compute  y := tau * A * v  storing y in TAU(1:i)
!
            CALL DSPMV( UPLO, I, TAUI, AP, AP( I1 ), 1, ZERO, TAU, 1)
!
!              Compute  w := y - 1/2 * tau * (y**T *v) * v
!
            ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, AP( I1 ), 1 )
            CALL DAXPY( I, ALPHA, AP( I1 ), 1, TAU, 1 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w**T - w * v**T
!
            CALL DSPR2( UPLO, I, -ONE, AP( I1 ), 1, TAU, 1, AP )
!
            AP( I1+I-1 ) = E( I )
         END IF
         D( I+1 ) = AP( I1+I )
         TAU( I ) = TAUI
         I1 = I1 - I
10    CONTINUE
      D( 1 ) = AP( 1 )
   ELSE
!
!        Reduce the lower triangle of A. II is the index in AP of
!        A(i,i) and I1I1 is the index of A(i+1,i+1).
!
      II = 1
      DO 20 I = 1, N - 1
         I1I1 = II + N - I + 1
!
!           Generate elementary reflector H(i) = I - tau * v * v**T
!           to annihilate A(i+2:n,i)
!
         CALL DLARFG( N-I, AP( II+1 ), AP( II+2 ), 1, TAUI )
         E( I ) = AP( II+1 )
!
         IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
            AP( II+1 ) = ONE
!
!              Compute  y := tau * A * v  storing y in TAU(i:n-1)
!
            CALL DSPMV( UPLO, N-I, TAUI, AP( I1I1 ), AP( II+1 ), 1,&
                 ZERO, TAU( I ), 1 )
!
!              Compute  w := y - 1/2 * tau * (y**T *v) * v
!
            ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1, AP( II+1 ), 1 )
            CALL DAXPY( N-I, ALPHA, AP( II+1 ), 1, TAU( I ), 1 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w**T - w * v**T
!
            CALL DSPR2( UPLO, N-I, -ONE, AP( II+1 ), 1, TAU( I ), 1,&
                 AP( I1I1 ) )
!
            AP( II+1 ) = E( I )
         END IF
         D( I ) = AP( II )
         TAU( I ) = TAUI
         II = I1I1
20    CONTINUE
      D( N ) = AP( II )
   END IF
!
   RETURN
!
!     End of DSPTRD
!
 END SUBROUTINE DSPTRD
!
!=
!
!> \brief \b DLARFG generates an elementary reflector (Householder matrix).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLARFG + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfg.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfg.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfg.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
! 
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       DOUBLE PRECISION   ALPHA, TAU
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   X( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARFG generates a real elementary reflector H of order n, such
!> that
!>
!>       H * ( alpha ) = ( beta ),   H**T * H = I.
!>           (   x   )   (   0  )
!>
!> where alpha and beta are scalars, and x is an (n-1)-element real
!> vector. H is represented in the form
!>
!>       H = I - tau * ( 1 ) * ( 1 v**T ) ,
!>                     ( v )
!>
!> where tau is a real scalar and v is a real (n-1)-element
!> vector.
!>
!> If the elements of x are all zero, then tau = 0 and H is taken to be
!> the unit matrix.
!>
!> Otherwise  1 <= tau <= 2.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the elementary reflector.
!> \endverbatim
!>
!> \param[in,out] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION
!>          On entry, the value alpha.
!>          On exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension
!>                         (1+(N-2)*abs(INCX))
!>          On entry, the vector x.
!>          On exit, it is overwritten with the vector v.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between elements of X. INCX > 0.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION
!>          The value tau.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup doubleOTHERauxiliary
!
!  =====================================================================
 SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
   INTEGER            INCX, N
   DOUBLE PRECISION   ALPHA, TAU
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE, ZERO
   PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
   INTEGER            J, KNT
   DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Functions ..
!?   DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
!   DOUBLE PRECISION   DNRM2
!   EXTERNAL           DNRM2
!     ..
!     .. Intrinsic Functions ..
!   INTRINSIC          ABS, SIGN
!     ..
!     .. External Subroutines ..
!   EXTERNAL           DSCAL
!     ..
!     .. Executable Statements ..
!
   IF( N.LE.1 ) THEN
      TAU = ZERO
      RETURN
   END IF
!
   XNORM = DNRM2( N-1, X, INCX )
!
   IF( XNORM.EQ.ZERO ) THEN
!
!        H  =  I
!
      TAU = ZERO
   ELSE
!
!        general case
!
      BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
      SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
      KNT = 0
      IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
         RSAFMN = ONE / SAFMIN
10       CONTINUE
         KNT = KNT + 1
         CALL DSCAL( N-1, RSAFMN, X, INCX )
         BETA = BETA*RSAFMN
         ALPHA = ALPHA*RSAFMN
         IF( ABS( BETA ).LT.SAFMIN )&
              GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
         XNORM = DNRM2( N-1, X, INCX )
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
      END IF
      TAU = ( BETA-ALPHA ) / BETA
      CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
!
!        If ALPHA is subnormal, it may lose relative accuracy
!
      DO 20 J = 1, KNT
         BETA = BETA*SAFMIN
20    CONTINUE
         ALPHA = BETA
   END IF
!
   RETURN
!
!     End of DLARFG
!
 END SUBROUTINE DLARFG
!
!=
!
!> \brief \b DLAPY2 returns sqrt(x2+y2).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAPY2 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy2.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy2.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy2.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   X, Y
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!> overflow.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is DOUBLE PRECISION
!>          X and Y specify the values x and y.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
 DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
   DOUBLE PRECISION   X, Y
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO
   PARAMETER          ( ZERO = 0.0D0 )
   DOUBLE PRECISION   ONE
   PARAMETER          ( ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
   DOUBLE PRECISION   W, XABS, YABS, Z
!     ..
!     .. Intrinsic Functions ..
!   INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
   XABS = ABS( X )
   YABS = ABS( Y )
   W = MAX( XABS, YABS )
   Z = MIN( XABS, YABS )
   IF( Z.EQ.ZERO ) THEN
      DLAPY2 = W
   ELSE
      DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
   END IF
   RETURN
!
!     End of DLAPY2
!
 END FUNCTION DLAPY2
!
!=
!
!> \brief \b DSTERF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DSTERF + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsterf.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsterf.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsterf.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSTERF( N, D, E, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), E( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
!> using the Pal-Walker-Kahan variant of the QL or QR algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          On entry, the n diagonal elements of the tridiagonal matrix.
!>          On exit, if INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          On entry, the (n-1) subdiagonal elements of the tridiagonal
!>          matrix.
!>          On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  the algorithm failed to find all of the eigenvalues in
!>                a total of 30*N iterations; if INFO = i, then i
!>                elements of E have not converged to zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
 SUBROUTINE DSTERF( N, D, E, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
   INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   D( * ), E( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE, TWO, THREE
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,&
        THREE = 3.0D0 )
   INTEGER            MAXIT
   PARAMETER          ( MAXIT = 30 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M,&
        NMAXIT
   DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC,&
        OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN,&
        SIGMA, SSFMAX, SSFMIN, RMAX
!     ..
!     .. External Functions ..
!   DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
!?      EXTERNAL           DLAMCH, DLANST, DLAPY2
!     ..
!     .. External Subroutines ..
!?   EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
!   EXTERNAL           DLAE2
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
!
!     Quick return if possible
!
   IF( N.LT.0 ) THEN
      INFO = -1
      CALL XERBLA( 'DSTERF', -INFO )
      RETURN
   END IF
   IF( N.LE.1 ) RETURN
!
!     Determine the unit roundoff for this environment.
!
   EPS = DLAMCH( 'E' )
   EPS2 = EPS**2
   SAFMIN = DLAMCH( 'S' )
   SAFMAX = ONE / SAFMIN
   SSFMAX = SQRT( SAFMAX ) / THREE
   SSFMIN = SQRT( SAFMIN ) / EPS2
   RMAX = DLAMCH( 'O' )
!
!     Compute the eigenvalues of the tridiagonal matrix.
!
   NMAXIT = N*MAXIT
   SIGMA = ZERO
   JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
   L1 = 1
!
10 CONTINUE
   IF( L1.GT.N ) GO TO 170
   IF( L1.GT.1 ) E( L1-1 ) = ZERO
   DO 20 M = L1, N - 1
      IF( ABS( E( M ) ).LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+&
           1 ) ) ) )*EPS ) THEN
         E( M ) = ZERO
         GO TO 30
      END IF
20 CONTINUE
   M = N
!
30 CONTINUE
   L = L1
   LSV = L
   LEND = M
   LENDSV = LEND
   L1 = M + 1
   IF( LEND.EQ.L ) GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
   ANORM = DLANST( 'M', LEND-L+1, D( L ), E( L ) )
   ISCALE = 0
   IF( ANORM.EQ.ZERO ) GO TO 10      
   IF( (ANORM.GT.SSFMAX) ) THEN
      ISCALE = 1
      CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, INFO )
      CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, INFO)
   ELSE IF( ANORM.LT.SSFMIN ) THEN
      ISCALE = 2
      CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, INFO)
      CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, INFO)
   END IF
!
   DO 40 I = L, LEND - 1
      E( I ) = E( I )**2
40 CONTINUE
!
!     Choose between QL and QR iteration
!
   IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
      LEND = LSV
      L = LENDSV
   END IF
!
   IF( LEND.GE.L ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
50    CONTINUE
      IF( L.NE.LEND ) THEN
         DO 60 M = L, LEND - 1
            IF( ABS( E( M ) ).LE.EPS2*ABS( D( M )*D( M+1 ) ) ) GO TO 70
60       CONTINUE
      END IF
      M = LEND
!
70    CONTINUE
      IF( M.LT.LEND ) E( M ) = ZERO
      P = D( L )
      IF( M.EQ.L ) GO TO 90
!
!        If remaining matrix is 2 by 2, use DLAE2 to compute its
!        eigenvalues.
!
      IF( M.EQ.L+1 ) THEN
         RTE = SQRT( E( L ) )
         CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
         D( L ) = RT1
         D( L+1 ) = RT2
         E( L ) = ZERO
         L = L + 2
         IF( L.LE.LEND ) GO TO 50
         GO TO 150
      END IF
!
      IF( JTOT.EQ.NMAXIT ) GO TO 150
      JTOT = JTOT + 1
!
!        Form shift.
!
      RTE = SQRT( E( L ) )
      SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
      R = DLAPY2( SIGMA, ONE )
      SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
!
      C = ONE
      S = ZERO
      GAMMA = D( M ) - SIGMA
      P = GAMMA*GAMMA
!
!        Inner loop
!
      DO 80 I = M - 1, L, -1
         BB = E( I )
         R = P + BB
         IF( I.NE.M-1 ) E( I+1 ) = S*R
         OLDC = C
         C = P / R
         S = BB / R
         OLDGAM = GAMMA
         ALPHA = D( I )
         GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
         D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
         IF( C.NE.ZERO ) THEN
            P = ( GAMMA*GAMMA ) / C
         ELSE
            P = OLDC*BB
         END IF
80    CONTINUE
!
      E( L ) = S*P
      D( L ) = SIGMA + GAMMA
      GO TO 50
!
!        Eigenvalue found.
!
90    CONTINUE
      D( L ) = P
!
      L = L + 1
      IF( L.LE.LEND ) GO TO 50
      GO TO 150
!
   ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
100   CONTINUE
      DO 110 M = L, LEND + 1, -1
         IF( ABS( E( M-1 ) ).LE.EPS2*ABS( D( M )*D( M-1 ) ) ) GO TO 120
110   CONTINUE
      M = LEND
!
120   CONTINUE
      IF( M.GT.LEND ) E( M-1 ) = ZERO
      P = D( L )
      IF( M.EQ.L ) GO TO 140
!
!        If remaining matrix is 2 by 2, use DLAE2 to compute its
!        eigenvalues.
!
      IF( M.EQ.L-1 ) THEN
         RTE = SQRT( E( L-1 ) )
         CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
         D( L ) = RT1
         D( L-1 ) = RT2
         E( L-1 ) = ZERO
         L = L - 2
         IF( L.GE.LEND ) GO TO 100
         GO TO 150
      END IF
!
      IF( JTOT.EQ.NMAXIT ) GO TO 150
      JTOT = JTOT + 1
!
!        Form shift.
!
      RTE = SQRT( E( L-1 ) )
      SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
      R = DLAPY2( SIGMA, ONE )
      SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
!
      C = ONE
      S = ZERO
      GAMMA = D( M ) - SIGMA
      P = GAMMA*GAMMA
!
!        Inner loop
!
      DO 130 I = M, L - 1
         BB = E( I )
         R = P + BB
         IF( I.NE.M ) E( I-1 ) = S*R
         OLDC = C
         C = P / R
         S = BB / R
         OLDGAM = GAMMA
         ALPHA = D( I+1 )
         GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
         D( I ) = OLDGAM + ( ALPHA-GAMMA )
         IF( C.NE.ZERO ) THEN
            P = ( GAMMA*GAMMA ) / C
         ELSE
            P = OLDC*BB
         END IF
130   CONTINUE
!
      E( L-1 ) = S*P
      D( L ) = SIGMA + GAMMA
      GO TO 100
!
!        Eigenvalue found.
!
140   CONTINUE
      D( L ) = P
!
      L = L - 1
      IF( L.GE.LEND ) GO TO 100
      GO TO 150
!
   END IF
!
!     Undo scaling if necessary
!
150 CONTINUE
   IF( ISCALE.EQ.1 ) &
        CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,&
        D( LSV ), N, INFO )
   IF( ISCALE.EQ.2 ) &
        CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,&
        D( LSV ), N, INFO )
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
   IF( JTOT.LT.NMAXIT ) GO TO 10
   DO 160 I = 1, N - 1
      IF( E( I ).NE.ZERO ) INFO = INFO + 1
160   CONTINUE
   GO TO 180
!
!     Sort eigenvalues in increasing order.
!
170 CONTINUE
   CALL DLASRT( 'I', N, D, INFO )
!
180 CONTINUE
   RETURN
!
!     End of DSTERF
!
 END SUBROUTINE DSTERF
!
!=
!
!> \brief \b DOPMTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DOPMTR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dopmtr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dopmtr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dopmtr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DOPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK,
!                          INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS, UPLO
!       INTEGER            INFO, LDC, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AP( * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DOPMTR overwrites the general real M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'T':      Q**T * C       C * Q**T
!>
!> where Q is a real orthogonal matrix of order nq, with nq = m if
!> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!> nq-1 elementary reflectors, as returned by DSPTRD using packed
!> storage:
!>
!> if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
!>
!> if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**T from the Left;
!>          = 'R': apply Q or Q**T from the Right.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U': Upper triangular packed storage used in previous
!>                 call to DSPTRD;
!>          = 'L': Lower triangular packed storage used in previous
!>                 call to DSPTRD.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'T':  Transpose, apply Q**T.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension
!>                               (M*(M+1)/2) if SIDE = 'L'
!>                               (N*(N+1)/2) if SIDE = 'R'
!>          The vectors which define the elementary reflectors, as
!>          returned by DSPTRD.  AP is modified by the routine but
!>          restored on exit.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (M-1) if SIDE = 'L'
!>                                     or (N-1) if SIDE = 'R'
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DSPTRD.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension
!>                                   (N) if SIDE = 'L'
!>                                   (M) if SIDE = 'R'
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
 SUBROUTINE DOPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK,&
      INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
   CHARACTER          SIDE, TRANS, UPLO
   INTEGER            INFO, LDC, M, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   AP( * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE
   PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            FORWRD, LEFT, NOTRAN, UPPER
   INTEGER            I, I1, I2, I3, IC, II, JC, MI, NI, NQ
   DOUBLE PRECISION   AII
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines .. DLARFG exists ...
!?      EXTERNAL           DLARF, XERBLA
!   EXTERNAL           DLARF
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
   INFO = 0
   LEFT = LSAME( SIDE, 'L' )
   NOTRAN = LSAME( TRANS, 'N' )
   UPPER = LSAME( UPLO, 'U' )
!
!     NQ is the order of Q
!
   IF( LEFT ) THEN
      NQ = M
   ELSE
      NQ = N
   END IF
   IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
      INFO = -1
   ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -2
   ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
      INFO = -3
   ELSE IF( M.LT.0 ) THEN
      INFO = -4
   ELSE IF( N.LT.0 ) THEN
      INFO = -5
   ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
      INFO = -9
   END IF
   IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'DOPMTR', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
!
   IF( UPPER ) THEN
!
!        Q was determined by a call to DSPTRD with UPLO = 'U'
!
      FORWRD = ( LEFT .AND. NOTRAN ) .OR.&
           ( .NOT.LEFT .AND. .NOT.NOTRAN )
!
      IF( FORWRD ) THEN
         I1 = 1
         I2 = NQ - 1
         I3 = 1
         II = 2
      ELSE
         I1 = NQ - 1
         I2 = 1
         I3 = -1
         II = NQ*( NQ+1 ) / 2 - 1
      END IF
!
      IF( LEFT ) THEN
         NI = N
      ELSE
         MI = M
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!              H(i) is applied to C(1:i,1:n)
!
            MI = I
         ELSE
!
!              H(i) is applied to C(1:m,1:i)
!
            NI = I
         END IF
!
!           Apply H(i)
!
         AII = AP( II )
         AP( II ) = ONE
         CALL DLARF( SIDE, MI, NI, AP( II-I+1 ), 1, TAU( I ), C, LDC, WORK )
         AP( II ) = AII
!
         IF( FORWRD ) THEN
            II = II + I + 2
         ELSE
            II = II - I - 1
         END IF
10    CONTINUE
   ELSE
!
!        Q was determined by a call to DSPTRD with UPLO = 'L'.
!
      FORWRD = ( LEFT .AND. .NOT.NOTRAN ) .OR.&
           ( .NOT.LEFT .AND. NOTRAN )
!
      IF( FORWRD ) THEN
         I1 = 1
         I2 = NQ - 1
         I3 = 1
         II = 2
      ELSE
         I1 = NQ - 1
         I2 = 1
         I3 = -1
         II = NQ*( NQ+1 ) / 2 - 1
      END IF
!
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
!
      DO 20 I = I1, I2, I3
         AII = AP( II )
         AP( II ) = ONE
         IF( LEFT ) THEN
!
!              H(i) is applied to C(i+1:m,1:n)
!
            MI = M - I
            IC = I + 1
         ELSE
!
!              H(i) is applied to C(1:m,i+1:n)
!
            NI = N - I
            JC = I + 1
         END IF
!
!           Apply H(i)
!
         CALL DLARF( SIDE, MI, NI, AP( II ), 1, TAU( I ),&
              C( IC, JC ), LDC, WORK )
         AP( II ) = AII
!
         IF( FORWRD ) THEN
            II = II + NQ - I + 1
         ELSE
            II = II - NQ + I - 2
         END IF
20    CONTINUE
   END IF
   RETURN
!
!     End of DOPMTR
!
 END SUBROUTINE DOPMTR
!
!=
!
!> \brief \b DLARF applies an elementary reflector to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLARF + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
! 
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            INCV, LDC, M, N
!       DOUBLE PRECISION   TAU
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARF applies a real elementary reflector H to a real m by n matrix
!> C, from either the left or the right. H is represented in the form
!>
!>       H = I - tau * v * v**T
!>
!> where tau is a real scalar and v is a real vector.
!>
!> If tau = 0, then H is taken to be the unit matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': form  H * C
!>          = 'R': form  C * H
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension
!>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!>          The vector v in the representation of H. V is not used if
!>          TAU = 0.
!> \endverbatim
!>
!> \param[in] INCV
!> \verbatim
!>          INCV is INTEGER
!>          The increment between elements of v. INCV <> 0.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
!>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!>          or C * H if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension
!>                         (N) if SIDE = 'L'
!>                      or (M) if SIDE = 'R'
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup doubleOTHERauxiliary
!
!  =====================================================================
 SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
   CHARACTER          SIDE
   INTEGER            INCV, LDC, M, N
   DOUBLE PRECISION   TAU
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE, ZERO
   PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            APPLYLEFT
   INTEGER            I, LASTV, LASTC
!     ..
!     .. External Subroutines ..
!   EXTERNAL           DGEMV, DGER
!   EXTERNAL           DGER
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!   INTEGER            ILADLR, ILADLC
!      EXTERNAL           LSAME, ILADLR, ILADLC
!   EXTERNAL           ILADLR, ILADLC
!     ..
!     .. Executable Statements ..
!
   APPLYLEFT = LSAME( SIDE, 'L' )
   LASTV = 0
   LASTC = 0
   IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
      IF( APPLYLEFT ) THEN
         LASTV = M
      ELSE
         LASTV = N
      END IF
      IF( INCV.GT.0 ) THEN
         I = 1 + (LASTV-1) * INCV
      ELSE
         I = 1
      END IF
!     Look for the last non-zero row in V.
      DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
         LASTV = LASTV - 1
         I = I - INCV
      END DO
      IF( APPLYLEFT ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
         LASTC = ILADLC(LASTV, N, C, LDC)
      ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
         LASTC = ILADLR(M, LASTV, C, LDC)
      END IF
   END IF
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
   IF( APPLYLEFT ) THEN
!
!        Form  H * C
!
      IF( LASTV.GT.0 ) THEN
!
!           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
!
         CALL DGEMV( 'Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV,&
              ZERO, WORK, 1 )
!
!           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
!
         CALL DGER( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
      END IF
   ELSE
!
!        Form  C * H
!
      IF( LASTV.GT.0 ) THEN
!
!           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
!
         CALL DGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC,&
              V, INCV, ZERO, WORK, 1 )
!
!           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
!
         CALL DGER( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
      END IF
   END IF
   RETURN
!
!     End of DLARF
!
 END SUBROUTINE DLARF
!
!=
!
!> \brief <b> DSPEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DSPEV + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dspev.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dspev.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dspev.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, LDZ, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSPEV computes all the eigenvalues and, optionally, eigenvectors of a
!> real symmetric matrix A in packed storage.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBZ
!> \verbatim
!>          JOBZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only;
!>          = 'V':  Compute eigenvalues and eigenvectors.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the symmetric matrix
!>          A, packed columnwise in a linear array.  The j-th column of A
!>          is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!>
!>          On exit, AP is overwritten by values generated during the
!>          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
!>          and first superdiagonal of the tridiagonal matrix T overwrite
!>          the corresponding elements of A, and if UPLO = 'L', the
!>          diagonal and first subdiagonal of T overwrite the
!>          corresponding elements of A.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
!>          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
!>          eigenvectors of the matrix A, with the i-th column of Z
!>          holding the eigenvector associated with W(i).
!>          If JOBZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          JOBZ = 'V', LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (3*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = i, the algorithm failed to converge; i
!>                off-diagonal elements of an intermediate tridiagonal
!>                form did not converge to zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup doubleOTHEReigen
!
!  =====================================================================
 SUBROUTINE DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
!
!  -- LAPACK driver routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
   CHARACTER          JOBZ, UPLO
   INTEGER            INFO, LDZ, N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            WANTZ
   INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE
   DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,&
        SMLNUM
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      DOUBLE PRECISION   DLAMCH, DLANSP
!      EXTERNAL           LSAME, DLAMCH, DLANSP
!     ..
!     .. External Subroutines ..
!   EXTERNAL           DOPGTR, DSCAL, DSPTRD, DSTEQR, DSTERF, XERBLA
!   EXTERNAL           DOPGTR, DSTEQR
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   WANTZ = LSAME( JOBZ, 'V' )
!
   INFO = 0
   IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
      INFO = -1
   ELSE IF( .NOT.( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) ) THEN
      INFO = -2
   ELSE IF( N.LT.0 ) THEN
      INFO = -3
   ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
      INFO = -7
   END IF
!
   IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'DSPEV ', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N.EQ.0 ) RETURN
!
   IF( N.EQ.1 ) THEN
      W( 1 ) = AP( 1 )
      IF( WANTZ ) Z( 1, 1 ) = ONE
      RETURN
   END IF
!
!     Get machine constants.
!
   SAFMIN = DLAMCH( 'Safe minimum' )
   EPS = DLAMCH( 'Precision' )
   SMLNUM = SAFMIN / EPS
   BIGNUM = ONE / SMLNUM
   RMIN = SQRT( SMLNUM )
   RMAX = SQRT( BIGNUM )
!
!     Scale matrix to allowable range, if necessary.
!
   ANRM = DLANSP( 'M', UPLO, N, AP, WORK )
   ISCALE = 0
   IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
      ISCALE = 1
      SIGMA = RMIN / ANRM
   ELSE IF( ANRM.GT.RMAX ) THEN
      ISCALE = 1
      SIGMA = RMAX / ANRM
   END IF
   IF( ISCALE.EQ.1 ) THEN
      CALL DSCAL( ( N*( N+1 ) ) / 2, SIGMA, AP, 1 )
   END IF
!
!     Call DSPTRD to reduce symmetric packed matrix to tridiagonal form.
!
   INDE = 1
   INDTAU = INDE + N
   CALL DSPTRD( UPLO, N, AP, W, WORK( INDE ), WORK( INDTAU ), IINFO )
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, first call
!     DOPGTR to generate the orthogonal matrix, then call DSTEQR.
!
   IF( .NOT.WANTZ ) THEN
      CALL DSTERF( N, W, WORK( INDE ), INFO )
   ELSE
      INDWRK = INDTAU + N
      CALL DOPGTR( UPLO, N, AP, WORK( INDTAU ), Z, LDZ,&
           WORK( INDWRK ), IINFO )
      CALL DSTEQR( JOBZ, N, W, WORK( INDE ), Z, LDZ, WORK( INDTAU ),&
           INFO )
   END IF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
   IF( ISCALE.EQ.1 ) THEN
      IF( INFO.EQ.0 ) THEN
         IMAX = N
      ELSE
         IMAX = INFO - 1
      END IF
      CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
   END IF
!
   RETURN
!
!     End of DSPEV
!
 END SUBROUTINE DSPEV
!
!=
!
!> \brief \b DSTEDC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DSTEDC + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dstedc.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dstedc.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dstedc.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK,
!                          LIWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          COMPZ
!       INTEGER            INFO, LDZ, LIWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSTEDC computes all eigenvalues and, optionally, eigenvectors of a
!> symmetric tridiagonal matrix using the divide and conquer method.
!> The eigenvectors of a full or band real symmetric matrix can also be
!> found if DSYTRD or DSPTRD or DSBTRD has been used to reduce this
!> matrix to tridiagonal form.
!>
!> This code makes very mild assumptions about floating point
!> arithmetic. It will work on machines with a guard digit in
!> add/subtract, or on those binary machines without guard digits
!> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!> It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.  See DLAED3 for details.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only.
!>          = 'I':  Compute eigenvectors of tridiagonal matrix also.
!>          = 'V':  Compute eigenvectors of original dense symmetric
!>                  matrix also.  On entry, Z contains the orthogonal
!>                  matrix used to reduce the original matrix to
!>                  tridiagonal form.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          On entry, the diagonal elements of the tridiagonal matrix.
!>          On exit, if INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          On entry, the subdiagonal elements of the tridiagonal matrix.
!>          On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ,N)
!>          On entry, if COMPZ = 'V', then Z contains the orthogonal
!>          matrix used in the reduction to tridiagonal form.
!>          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
!>          orthonormal eigenvectors of the original symmetric matrix,
!>          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!>          of the symmetric tridiagonal matrix.
!>          If  COMPZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1.
!>          If eigenvectors are desired, then LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array,
!>                                         dimension (LWORK)
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If COMPZ = 'N' or N <= 1 then LWORK must be at least 1.
!>          If COMPZ = 'V' and N > 1 then LWORK must be at least
!>                         ( 1 + 3*N + 2*N*lg N + 4*N**2 ),
!>                         where lg( N ) = smallest integer k such
!>                         that 2**k >= N.
!>          If COMPZ = 'I' and N > 1 then LWORK must be at least
!>                         ( 1 + 4*N + N**2 ).
!>          Note that for COMPZ = 'I' or 'V', then if N is less than or
!>          equal to the minimum divide size, usually 25, then LWORK need
!>          only be max(1,2*(N-1)).
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1.
!>          If COMPZ = 'V' and N > 1 then LIWORK must be at least
!>                         ( 6 + 6*N + 5*N*lg N ).
!>          If COMPZ = 'I' and N > 1 then LIWORK must be at least
!>                         ( 3 + 5*N ).
!>          Note that for COMPZ = 'I' or 'V', then if N is less than or
!>          equal to the minimum divide size, usually 25, then LIWORK
!>          need only be 1.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of the IWORK array,
!>          returns this value as the first entry of the IWORK array, and
!>          no error message related to LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  The algorithm failed to compute an eigenvalue while
!>                working on the submatrix lying in rows and columns
!>                INFO/(N+1) through mod(INFO,N+1).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2015
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA \n
!>  Modified by Francoise Tisseur, University of Tennessee
!>
!  =====================================================================
 SUBROUTINE DSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK,&
      LIWORK, INFO )
!
!  -- LAPACK computational routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
   CHARACTER          COMPZ
   INTEGER            INFO, LDZ, LIWORK, LWORK, N
!     ..
!     .. Array Arguments ..
   INTEGER            IWORK( * )
   DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE, TWO
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            LQUERY
   INTEGER            FINISH, I, ICOMPZ, II, J, K, LGN, LIWMIN,&
        LWMIN, M, SMLSIZ, START, STOREZ, STRTRW
   DOUBLE PRECISION   EPS, ORGNRM, P, TINY
!     ..
!     .. External Functions ..
!   LOGICAL            LSAME
!   INTEGER            ILAENV
!   DOUBLE PRECISION  DLANST
!      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANST
!   EXTERNAL           DLANST
!     ..
!     .. External Subroutines ..
!   EXTERNAL           DGEMM, DLACPY, DLAED0, DLASCL, DLASET, DLASRT,&
!        DSTEQR, DSTERF, DSWAP, XERBLA
!   EXTERNAL           DSTEQR
!     ..
!     .. Intrinsic Functions ..
!   INTRINSIC          ABS, DBLE, INT, LOG, MAX, MOD, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
   LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
!
   IF( LSAME( COMPZ, 'N' ) ) THEN
      ICOMPZ = 0
   ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
      ICOMPZ = 1
   ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
      ICOMPZ = 2
   ELSE
      ICOMPZ = -1
   END IF
   IF( ICOMPZ.LT.0 ) THEN
      INFO = -1
   ELSE IF( N.LT.0 ) THEN
      INFO = -2
   ELSE IF( ( LDZ.LT.1 ) .OR.&
        ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
      INFO = -6
   END IF
!
   IF( INFO.EQ.0 ) THEN
!
!        Compute the workspace requirements
!
      SMLSIZ = ILAENV( 9, 'DSTEDC', ' ', 0, 0, 0, 0 )
      IF( N.LE.1 .OR. ICOMPZ.EQ.0 ) THEN
         LIWMIN = 1
         LWMIN = 1
      ELSE IF( N.LE.SMLSIZ ) THEN
         LIWMIN = 1
         LWMIN = 2*( N - 1 )
      ELSE
         LGN = INT( LOG( DBLE( N ) )/LOG( TWO ) )
         IF( 2**LGN.LT.N ) LGN = LGN + 1
         IF( 2**LGN.LT.N ) LGN = LGN + 1
         IF( ICOMPZ.EQ.1 ) THEN
            LWMIN = 1 + 3*N + 2*N*LGN + 4*N**2
            LIWMIN = 6 + 6*N + 5*N*LGN
         ELSE IF( ICOMPZ.EQ.2 ) THEN
            LWMIN = 1 + 4*N + N**2
            LIWMIN = 3 + 5*N
         END IF
      END IF
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
!
      IF( LWORK.LT.LWMIN .AND. .NOT. LQUERY ) THEN
         INFO = -8
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT. LQUERY ) THEN
         INFO = -10
      END IF
   END IF
!
   IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'DSTEDC', -INFO )
      RETURN
   ELSE IF (LQUERY) THEN
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N.EQ.0 ) RETURN
   IF( N.EQ.1 ) THEN
      IF( ICOMPZ.NE.0 ) Z( 1, 1 ) = ONE
      RETURN
   END IF
!
!     If the following conditional clause is removed, then the routine
!     will use the Divide and Conquer routine to compute only the
!     eigenvalues, which requires (3N + 3N**2) real workspace and
!     (2 + 5N + 2N lg(N)) integer workspace.
!     Since on many architectures DSTERF is much faster than any other
!     algorithm for finding eigenvalues only, it is used here
!     as the default. If the conditional clause is removed, then
!     information on the size of workspace needs to be changed.
!
!     If COMPZ = 'N', use DSTERF to compute the eigenvalues.
!
   IF( ICOMPZ.EQ.0 ) THEN
      CALL DSTERF( N, D, E, INFO )
      GO TO 50
   END IF
!
!     If N is smaller than the minimum divide size (SMLSIZ+1), then
!     solve the problem with another solver.
!
   IF( N.LE.SMLSIZ ) THEN
!
      CALL DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
   ELSE
!
!        If COMPZ = 'V', the Z matrix must be stored elsewhere for later
!        use.
!
      IF( ICOMPZ.EQ.1 ) THEN
         STOREZ = 1 + N*N
      ELSE
         STOREZ = 1
      END IF
!
      IF( ICOMPZ.EQ.2 ) THEN
         CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
      END IF
!
!        Scale.
!
      ORGNRM = DLANST( 'M', N, D, E )
      IF( ORGNRM.EQ.ZERO ) GO TO 50
!
      EPS = DLAMCH( 'Epsilon' )
!
      START = 1
!
!        while ( START <= N )
!
10    CONTINUE
      IF( START.LE.N ) THEN
!
!           Let FINISH be the position of the next subdiagonal entry
!           such that E( FINISH ) <= TINY or FINISH = N if no such
!           subdiagonal exists.  The matrix identified by the elements
!           between START and FINISH constitutes an independent
!           sub-problem.
!
         FINISH = START
20       CONTINUE
         IF( FINISH.LT.N ) THEN
            TINY = EPS*SQRT( ABS( D( FINISH ) ) )*&
                 SQRT( ABS( D( FINISH+1 ) ) )
            IF( ABS( E( FINISH ) ).GT.TINY ) THEN
               FINISH = FINISH + 1
               GO TO 20
            END IF
         END IF
!
!           (Sub) Problem determined.  Compute its size and solve it.
!
         M = FINISH - START + 1
         IF( M.EQ.1 ) THEN
            START = FINISH + 1
            GO TO 10
         END IF
         IF( M.GT.SMLSIZ ) THEN
!
!              Scale.
!
            ORGNRM = DLANST( 'M', M, D( START ), E( START ) )
            CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M, 1, D( START ), M,INFO)
            CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M-1, 1, E( START ),M-1, INFO )
!
            IF( ICOMPZ.EQ.1 ) THEN
               STRTRW = 1
            ELSE
               STRTRW = START
            END IF
            CALL DLAED0( ICOMPZ, N, M, D( START ), E( START ),&
                 Z( STRTRW, START ), LDZ, WORK( 1 ), N,&
                 WORK( STOREZ ), IWORK, INFO )
            IF( INFO.NE.0 ) THEN
               INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) +&
                    MOD( INFO, ( M+1 ) ) + START - 1
               GO TO 50
            END IF
!
!              Scale back.
!
            CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, M, 1, D( START ), M,INFO)
!
         ELSE
            IF( ICOMPZ.EQ.1 ) THEN
!
!                 Since QR won't update a Z matrix which is larger than
!                 the length of D, we must solve the sub-problem in a
!                 workspace and then multiply back into Z.
!
               CALL DSTEQR( 'I', M, D( START ), E( START ), WORK, M,&
                    WORK( M*M+1 ), INFO )
               CALL DLACPY( 'A', N, M, Z( 1, START ), LDZ,&
                    WORK( STOREZ ), N )
               CALL DGEMM( 'N', 'N', N, M, M, ONE,&
                    WORK( STOREZ ), N, WORK, M, ZERO,&
                    Z( 1, START ), LDZ )
            ELSE IF( ICOMPZ.EQ.2 ) THEN
               CALL DSTEQR( 'I', M, D( START ), E( START ),&
                    Z( START, START ), LDZ, WORK, INFO )
            ELSE
               CALL DSTERF( M, D( START ), E( START ), INFO )
            END IF
            IF( INFO.NE.0 ) THEN
               INFO = START*( N+1 ) + FINISH
               GO TO 50
            END IF
         END IF
!
         START = FINISH + 1
         GO TO 10
      END IF
!
!        endwhile
!
      IF( ICOMPZ.EQ.0 ) THEN
!
!          Use Quick Sort
!
         CALL DLASRT( 'I', N, D, INFO )
!
      ELSE
!
!          Use Selection Sort to minimize swaps of eigenvectors
!
         DO 40 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 30 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
30          CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
40       CONTINUE
         END IF
      END IF
!
50    CONTINUE
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
!
      RETURN
!
!     End of DSTEDC
!
   END SUBROUTINE DSTEDC
!
!=
!
!> \brief \b DLANST returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a real symmetric tridiagonal matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLANST + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanst.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanst.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanst.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
! 
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), E( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLANST  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> real symmetric tridiagonal matrix A.
!> \endverbatim
!>
!> \return DLANST
!> \verbatim
!>
!>    DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!>             (
!>             ( norm1(A),         NORM = '1', 'O' or 'o'
!>             (
!>             ( normI(A),         NORM = 'I' or 'i'
!>             (
!>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normF  denotes the  Frobenius norm of a matrix (square root of sum of
!> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies the value to be returned in DLANST as described
!>          above.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.  When N = 0, DLANST is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal elements of A.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) sub-diagonal or super-diagonal elements of A.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
   DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
     CHARACTER          NORM
     INTEGER            N
!     ..
!     .. Array Arguments ..
     DOUBLE PRECISION   D( * ), E( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
     DOUBLE PRECISION   ONE, ZERO
     PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
     INTEGER            I
     DOUBLE PRECISION   ANORM, SCALE, SUM
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME, DISNAN
!      EXTERNAL           LSAME, DISNAN
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLASSQ
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
     IF( N.LE.0 ) THEN
        ANORM = ZERO
     ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
        ANORM = ABS( D( N ) )
        DO 10 I = 1, N - 1
           SUM = ABS( D( I ) )
           IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
           SUM = ABS( E( I ) )
           IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
10      CONTINUE
     ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.&
          LSAME( NORM, 'I' ) ) THEN
!
!        Find norm1(A).
!
        IF( N.EQ.1 ) THEN
           ANORM = ABS( D( 1 ) )
        ELSE
           ANORM = ABS( D( 1 ) )+ABS( E( 1 ) )
           SUM = ABS( E( N-1 ) )+ABS( D( N ) )
           IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
           DO 20 I = 2, N - 1
              SUM = ABS( D( I ) )+ABS( E( I ) )+ABS( E( I-1 ) )
              IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
20         CONTINUE
        END IF
     ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
        SCALE = ZERO
        SUM = ONE
        IF( N.GT.1 ) THEN
           CALL DLASSQ( N-1, E, 1, SCALE, SUM )
           SUM = 2*SUM
        END IF
        CALL DLASSQ( N, D, 1, SCALE, SUM )
        ANORM = SCALE*SQRT( SUM )
     END IF
!
     DLANST = ANORM
     RETURN
!
!     End of DLANST
!
  END FUNCTION DLANST
!
!=
!
!> \brief \b DLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLASCL + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlascl.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlascl.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlascl.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          TYPE
!       INTEGER            INFO, KL, KU, LDA, M, N
!       DOUBLE PRECISION   CFROM, CTO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASCL multiplies the M by N real matrix A by the real scalar
!> CTO/CFROM.  This is done without over/underflow as long as the final
!> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!> A may be full, upper triangular, lower triangular, upper Hessenberg,
!> or banded.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TYPE
!> \verbatim
!>          TYPE is CHARACTER*1
!>          TYPE indices the storage type of the input matrix.
!>          = 'G':  A is a full matrix.
!>          = 'L':  A is a lower triangular matrix.
!>          = 'U':  A is an upper triangular matrix.
!>          = 'H':  A is an upper Hessenberg matrix.
!>          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the lower
!>                  half stored.
!>          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the upper
!>                  half stored.
!>          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!>                  bandwidth KU. See DGBTRF for storage details.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] CFROM
!> \verbatim
!>          CFROM is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] CTO
!> \verbatim
!>          CTO is DOUBLE PRECISION
!>
!>          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!>          without over/underflow if the final result CTO*A(I,J)/CFROM
!>          can be represented without over/underflow.  CFROM must be
!>          nonzero.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!>          storage type.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          0  - successful exit
!>          <0 - if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
  SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
    CHARACTER          TYPE
    INTEGER            INFO, KL, KU, LDA, M, N
    DOUBLE PRECISION   CFROM, CTO
!     ..
!     .. Array Arguments ..
    DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
    DOUBLE PRECISION   ZERO, ONE
    PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
    LOGICAL            DONE
    INTEGER            I, ITYPE, J, K1, K2, K3, K4
    DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME, DISNAN
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           LSAME, DLAMCH, DISNAN
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. External Subroutines ..
!    EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
    INFO = 0
!
    IF( LSAME( TYPE, 'G' ) ) THEN
       ITYPE = 0
    ELSE IF( LSAME( TYPE, 'L' ) ) THEN
       ITYPE = 1
    ELSE IF( LSAME( TYPE, 'U' ) ) THEN
       ITYPE = 2
    ELSE IF( LSAME( TYPE, 'H' ) ) THEN
       ITYPE = 3
    ELSE IF( LSAME( TYPE, 'B' ) ) THEN
       ITYPE = 4
    ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
       ITYPE = 5
    ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
       ITYPE = 6
    ELSE
       ITYPE = -1
    END IF
!
    IF( ITYPE.EQ.-1 ) THEN
       INFO = -1
    ELSE IF( CFROM.EQ.ZERO .OR. DISNAN(CFROM) ) THEN
       INFO = -4
    ELSE IF( DISNAN(CTO) ) THEN
       INFO = -5
    ELSE IF( M.LT.0 ) THEN
       INFO = -6
    ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.&
         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
       INFO = -7
    ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
       INFO = -9
    ELSE IF( ITYPE.GE.4 ) THEN
       IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
          INFO = -2
       ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.&
            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )&
            THEN
          INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.&
              ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.&
              ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. M.EQ.0 ) RETURN
!
!     Get machine parameters
!
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
!
      CFROMC = CFROM
      CTOC = CTO
!
10    CONTINUE
      CFROM1 = CFROMC*SMLNUM
      IF( CFROM1.EQ.CFROMC ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         MUL = CTOC / CFROMC
         DONE = .TRUE.
         CTO1 = CTOC
      ELSE
         CTO1 = CTOC / BIGNUM
         IF( CTO1.EQ.CTOC ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            MUL = CTOC
            DONE = .TRUE.
            CFROMC = ONE
         ELSE IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
            MUL = SMLNUM
            DONE = .FALSE.
            CFROMC = CFROM1
         ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
            MUL = BIGNUM
            DONE = .FALSE.
            CTOC = CTO1
         ELSE
            MUL = CTOC / CFROMC
            DONE = .TRUE.
         END IF
      END IF
!
      IF( ITYPE.EQ.0 ) THEN
!
!        Full matrix
!
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
!
      ELSE IF( ITYPE.EQ.1 ) THEN
!
!        Lower triangular matrix
!
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
!
      ELSE IF( ITYPE.EQ.2 ) THEN
!
!        Upper triangular matrix
!
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
!
      ELSE IF( ITYPE.EQ.3 ) THEN
!
!        Upper Hessenberg matrix
!
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
!
      ELSE IF( ITYPE.EQ.4 ) THEN
!
!        Lower half of a symmetric band matrix
!
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
!
      ELSE IF( ITYPE.EQ.5 ) THEN
!
!        Upper half of a symmetric band matrix
!
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
!
      ELSE IF( ITYPE.EQ.6 ) THEN
!
!        Band matrix
!
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
!
      END IF
!
      IF( .NOT.DONE ) GO TO 10
!
      RETURN
!
!     End of DLASCL
!
   END SUBROUTINE DLASCL
!
!=
!
!> \brief \b DSTEQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DSTEQR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsteqr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsteqr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsteqr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          COMPZ
!       INTEGER            INFO, LDZ, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!> symmetric tridiagonal matrix using the implicit QL or QR method.
!> The eigenvectors of a full or band symmetric matrix can also be found
!> if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
!> tridiagonal form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only.
!>          = 'V':  Compute eigenvalues and eigenvectors of the original
!>                  symmetric matrix.  On entry, Z must contain the
!>                  orthogonal matrix used to reduce the original matrix
!>                  to tridiagonal form.
!>          = 'I':  Compute eigenvalues and eigenvectors of the
!>                  tridiagonal matrix.  Z is initialized to the identity
!>                  matrix.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          On entry, the diagonal elements of the tridiagonal matrix.
!>          On exit, if INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          On entry, the (n-1) subdiagonal elements of the tridiagonal
!>          matrix.
!>          On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
!>          On entry, if  COMPZ = 'V', then Z contains the orthogonal
!>          matrix used in the reduction to tridiagonal form.
!>          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
!>          orthonormal eigenvectors of the original symmetric matrix,
!>          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!>          of the symmetric tridiagonal matrix.
!>          If COMPZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          eigenvectors are desired, then  LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (max(1,2*N-2))
!>          If COMPZ = 'N', then WORK is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  the algorithm has failed to find all the eigenvalues in
!>                a total of 30*N iterations; if INFO = i, then i
!>                elements of E have not converged to zero; on exit, D
!>                and E contain the elements of a symmetric tridiagonal
!>                matrix which is orthogonally similar to the original
!>                matrix.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
   SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
     CHARACTER          COMPZ
     INTEGER            INFO, LDZ, N
!     ..
!     .. Array Arguments ..
     DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
     DOUBLE PRECISION   ZERO, ONE, TWO, THREE
     PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,&
          THREE = 3.0D0 )
     INTEGER            MAXIT
     PARAMETER          ( MAXIT = 30 )
!     ..
!     .. Local Scalars ..
     INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,&
          LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,&
          NM1, NMAXIT
     DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,&
          S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
!     ..
!     .. External Functions ..
!     LOGICAL            LSAME
!     DOUBLE PRECISION   DLAPY2
!      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASET, DLASR,
!     $                   DLASRT, DSWAP, XERBLA
!     EXTERNAL           DLARTG
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
     INFO = 0
!
     IF( LSAME( COMPZ, 'N' ) ) THEN
        ICOMPZ = 0
     ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
        ICOMPZ = 1
     ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
        ICOMPZ = 2
     ELSE
        ICOMPZ = -1
     END IF
     IF( ICOMPZ.LT.0 ) THEN
        INFO = -1
     ELSE IF( N.LT.0 ) THEN
        INFO = -2
     ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,&
          N ) ) ) THEN
        INFO = -6
     END IF
     IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DSTEQR', -INFO )
        RETURN
     END IF
!
!     Quick return if possible
!
     IF( N.EQ.0 ) RETURN
!
     IF( N.EQ.1 ) THEN
        IF( ICOMPZ.EQ.2 ) Z( 1, 1 ) = ONE
        RETURN
     END IF
!
!     Determine the unit roundoff and over/underflow thresholds.
!
     EPS = DLAMCH( 'E' )
     EPS2 = EPS**2
     SAFMIN = DLAMCH( 'S' )
     SAFMAX = ONE / SAFMIN
     SSFMAX = SQRT( SAFMAX ) / THREE
     SSFMIN = SQRT( SAFMIN ) / EPS2
!
!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
     IF( ICOMPZ.EQ.2 )&
          CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
!
     NMAXIT = N*MAXIT
     JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
     L1 = 1
     NM1 = N - 1
!
10   CONTINUE
     IF( L1.GT.N ) GO TO 160
     IF( L1.GT.1 ) E( L1-1 ) = ZERO
     IF( L1.LE.NM1 ) THEN
        DO 20 M = L1, NM1
           TST = ABS( E( M ) )
           IF( TST.EQ.ZERO ) GO TO 30
           IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+&
                1 ) ) ) )*EPS ) THEN
              E( M ) = ZERO
              GO TO 30
           END IF
20      CONTINUE
     END IF
     M = N
!
30   CONTINUE
     L = L1
     LSV = L
     LEND = M
     LENDSV = LEND
     L1 = M + 1
     IF( LEND.EQ.L ) GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
     ANORM = DLANST( 'M', LEND-L+1, D( L ), E( L ) )
     ISCALE = 0
     IF( ANORM.EQ.ZERO ) GO TO 10
     IF( ANORM.GT.SSFMAX ) THEN
        ISCALE = 1
        CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, INFO )
        CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, INFO )
     ELSE IF( ANORM.LT.SSFMIN ) THEN
        ISCALE = 2
        CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, INFO )
        CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, INFO )
     END IF
!
!     Choose between QL and QR iteration
!
     IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
        LEND = LSV
        L = LENDSV
     END IF
!
     IF( LEND.GT.L ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
40      CONTINUE
        IF( L.NE.LEND ) THEN
           LENDM1 = LEND - 1
           DO 50 M = L, LENDM1
              TST = ABS( E( M ) )**2
              IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+SAFMIN )GO TO 60
50         CONTINUE
        END IF
!
        M = LEND
!
60      CONTINUE
        IF( M.LT.LEND ) E( M ) = ZERO
        P = D( L )
        IF( M.EQ.L ) GO TO 80
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
        IF( M.EQ.L+1 ) THEN
           IF( ICOMPZ.GT.0 ) THEN
              CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
              WORK( L ) = C
              WORK( N-1+L ) = S
              CALL DLASR( 'R', 'V', 'B', N, 2, WORK( L ),&
                   WORK( N-1+L ), Z( 1, L ), LDZ )
           ELSE
              CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
           END IF
           D( L ) = RT1
           D( L+1 ) = RT2
           E( L ) = ZERO
           L = L + 2
           IF( L.LE.LEND ) GO TO 40
           GO TO 140
        END IF
!
        IF( JTOT.EQ.NMAXIT ) GO TO 140
        JTOT = JTOT + 1
!
!        Form shift.
!
        G = ( D( L+1 )-P ) / ( TWO*E( L ) )
        R = DLAPY2( G, ONE )
        G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
!
        S = ONE
        C = ONE
        P = ZERO
!
!        Inner loop
!
        MM1 = M - 1
        DO 70 I = MM1, L, -1
           F = S*E( I )
           B = C*E( I )
           CALL DLARTG( G, F, C, S, R )
           IF( I.NE.M-1 ) E( I+1 ) = R
           G = D( I+1 ) - P
           R = ( D( I )-G )*S + TWO*C*B
           P = S*R
           D( I+1 ) = G + P
           G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
           IF( ICOMPZ.GT.0 ) THEN
              WORK( I ) = C
              WORK( N-1+I ) = -S
           END IF
!
70      CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
        IF( ICOMPZ.GT.0 ) THEN
           MM = M - L + 1
           CALL DLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),&
                Z( 1, L ), LDZ )
        END IF
!
        D( L ) = D( L ) - P
        E( L ) = G
        GO TO 40
!
!        Eigenvalue found.
!
80      CONTINUE
        D( L ) = P
!
        L = L + 1
        IF( L.LE.LEND ) GO TO 40
        GO TO 140
!
     ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
90      CONTINUE
        IF( L.NE.LEND ) THEN
           LENDP1 = LEND + 1
           DO 100 M = L, LENDP1, -1
              TST = ABS( E( M-1 ) )**2
              IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+&
                   SAFMIN )GO TO 110
100        CONTINUE
        END IF
!
        M = LEND
!
110     CONTINUE
        IF( M.GT.LEND ) E( M-1 ) = ZERO
        P = D( L )
        IF( M.EQ.L ) GO TO 130
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
        IF( M.EQ.L-1 ) THEN
           IF( ICOMPZ.GT.0 ) THEN
              CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
              WORK( M ) = C
              WORK( N-1+M ) = S
              CALL DLASR( 'R', 'V', 'F', N, 2, WORK( M ), &
                   WORK( N-1+M ), Z( 1, L-1 ), LDZ )
           ELSE
              CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
           END IF
           D( L-1 ) = RT1
           D( L ) = RT2
           E( L-1 ) = ZERO
           L = L - 2
           IF( L.GE.LEND ) GO TO 90
           GO TO 140
        END IF
!
        IF( JTOT.EQ.NMAXIT ) GO TO 140
        JTOT = JTOT + 1
!
!        Form shift.
!
        G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
        R = DLAPY2( G, ONE )
        G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
!
        S = ONE
        C = ONE
        P = ZERO
!
!        Inner loop
!
        LM1 = L - 1
        DO 120 I = M, LM1
           F = S*E( I )
           B = C*E( I )
           CALL DLARTG( G, F, C, S, R )
           IF( I.NE.M ) E( I-1 ) = R
           G = D( I ) - P
           R = ( D( I+1 )-G )*S + TWO*C*B
           P = S*R
           D( I ) = G + P
           G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
           IF( ICOMPZ.GT.0 ) THEN
              WORK( I ) = C
              WORK( N-1+I ) = S
           END IF
!
120     CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
        IF( ICOMPZ.GT.0 ) THEN
           MM = L - M + 1
           CALL DLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),&
                Z( 1, M ), LDZ )
        END IF
!
        D( L ) = D( L ) - P
        E( LM1 ) = G
        GO TO 90
!
!        Eigenvalue found.
!
130     CONTINUE
        D( L ) = P
!
        L = L - 1
        IF( L.GE.LEND ) GO TO 90
        GO TO 140
!
     END IF
!
!     Undo scaling if necessary
!
140  CONTINUE
     IF( ISCALE.EQ.1 ) THEN
        CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,&
             D( LSV ), N, INFO )
        CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),&
             N, INFO )
     ELSE IF( ISCALE.EQ.2 ) THEN
        CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,&
             D( LSV ), N, INFO )
        CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),&
             N, INFO )
     END IF
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
     IF( JTOT.LT.NMAXIT ) GO TO 10
     DO 150 I = 1, N - 1
        IF( E( I ).NE.ZERO ) INFO = INFO + 1
150  CONTINUE
     GO TO 190
!
!     Order eigenvalues and eigenvectors.
!
160  CONTINUE
     IF( ICOMPZ.EQ.0 ) THEN
!
!        Use Quick Sort
!
        CALL DLASRT( 'I', N, D, INFO )
!
     ELSE
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
        DO 180 II = 2, N
           I = II - 1
           K = I
           P = D( I )
           DO 170 J = II, N
              IF( D( J ).LT.P ) THEN
                 K = J
                 P = D( J )
              END IF
170        CONTINUE
           IF( K.NE.I ) THEN
              D( K ) = D( I )
              D( I ) = P
              CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
           END IF
180     CONTINUE
     END IF
!
190  CONTINUE
     RETURN
!
!     End of DSTEQR
!
   END SUBROUTINE DSTEQR
!
!=
!
!> \brief \b DLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given values.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLASET + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaset.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaset.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaset.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, M, N
!       DOUBLE PRECISION   ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASET initializes an m-by-n matrix A to BETA on the diagonal and
!> ALPHA on the offdiagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be set.
!>          = 'U':      Upper triangular part is set; the strictly lower
!>                      triangular part of A is not changed.
!>          = 'L':      Lower triangular part is set; the strictly upper
!>                      triangular part of A is not changed.
!>          Otherwise:  All of the matrix A is set.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION
!>          The constant to which the offdiagonal elements are to be set.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION
!>          The constant to which the diagonal elements are to be set.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On exit, the leading m-by-n submatrix of A is set as follows:
!>
!>          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!>          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!>          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!>
!>          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2015
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
   SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
     CHARACTER          UPLO
     INTEGER            LDA, M, N
     DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. Array Arguments ..
     DOUBLE PRECISION   A( LDA, * )
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
10          CONTINUE
20       CONTINUE
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!
!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
!
      ELSE
!
!        Set the leading m-by-n submatrix to ALPHA.
!
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
!
!     Set the first min(M,N) diagonal elements to BETA.
!
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
!
      RETURN
!
!     End of DLASET
!
    END SUBROUTINE DLASET
!
!=
!
!> \brief \b DLAEV2 computes the eigenvalues and eigenvectors of a 2-by-2 symmetric/Hermitian matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAEV2 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaev2.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaev2.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaev2.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
!>    [  A   B  ]
!>    [  B   C  ].
!> On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
!> eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
!> eigenvector for RT1, giving the decomposition
!>
!>    [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
!>    [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION
!>          The (1,1) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION
!>          The (1,2) element and the conjugate of the (2,1) element of
!>          the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!>          The (2,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[out] RT1
!> \verbatim
!>          RT1 is DOUBLE PRECISION
!>          The eigenvalue of larger absolute value.
!> \endverbatim
!>
!> \param[out] RT2
!> \verbatim
!>          RT2 is DOUBLE PRECISION
!>          The eigenvalue of smaller absolute value.
!> \endverbatim
!>
!> \param[out] CS1
!> \verbatim
!>          CS1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] SN1
!> \verbatim
!>          SN1 is DOUBLE PRECISION
!>          The vector (CS1, SN1) is a unit right eigenvector for RT1.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  RT1 is accurate to a few ulps barring over/underflow.
!>
!>  RT2 may be inaccurate if there is massive cancellation in the
!>  determinant A*C-B*B; higher precision or correctly rounded or
!>  correctly truncated arithmetic would be needed to compute RT2
!>  accurately in all cases.
!>
!>  CS1 and SN1 are accurate to a few ulps barring over/underflow.
!>
!>  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!>  Underflow is harmless if the input data is 0 or exceeds
!>     underflow_threshold / macheps.
!> \endverbatim
!>
!  =====================================================================
    SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,&
           TB, TN
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
!     Compute the eigenvalues
!
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
!
!        Includes case AB=ADF=0
!
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
!
!        Includes case RT1 = RT2 = 0
!
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
!
!     Compute the eigenvector
!
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
!
!     End of DLAEV2
!
    END SUBROUTINE DLAEV2
!
!> \brief \b DLASR applies a sequence of plane rotations to a general rectangular matrix.
!
!=
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLASR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
! 
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, PIVOT, SIDE
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASR applies a sequence of plane rotations to a real matrix A,
!> from either the left or the right.
!> 
!> When SIDE = 'L', the transformation takes the form
!> 
!>    A := P*A
!> 
!> and when SIDE = 'R', the transformation takes the form
!> 
!>    A := A*P**T
!> 
!> where P is an orthogonal matrix consisting of a sequence of z plane
!> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
!> and P**T is the transpose of P.
!> 
!> When DIRECT = 'F' (Forward sequence), then
!> 
!>    P = P(z-1) * ... * P(2) * P(1)
!> 
!> and when DIRECT = 'B' (Backward sequence), then
!> 
!>    P = P(1) * P(2) * ... * P(z-1)
!> 
!> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
!> 
!>    R(k) = (  c(k)  s(k) )
!>         = ( -s(k)  c(k) ).
!> 
!> When PIVOT = 'V' (Variable pivot), the rotation is performed
!> for the plane (k,k+1), i.e., P(k) has the form
!> 
!>    P(k) = (  1                                            )
!>           (       ...                                     )
!>           (              1                                )
!>           (                   c(k)  s(k)                  )
!>           (                  -s(k)  c(k)                  )
!>           (                                1              )
!>           (                                     ...       )
!>           (                                            1  )
!> 
!> where R(k) appears as a rank-2 modification to the identity matrix in
!> rows and columns k and k+1.
!> 
!> When PIVOT = 'T' (Top pivot), the rotation is performed for the
!> plane (1,k+1), so P(k) has the form
!> 
!>    P(k) = (  c(k)                    s(k)                 )
!>           (         1                                     )
!>           (              ...                              )
!>           (                     1                         )
!>           ( -s(k)                    c(k)                 )
!>           (                                 1             )
!>           (                                      ...      )
!>           (                                             1 )
!> 
!> where R(k) appears in rows and columns 1 and k+1.
!> 
!> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
!> performed for the plane (k,z), giving P(k) the form
!> 
!>    P(k) = ( 1                                             )
!>           (      ...                                      )
!>           (             1                                 )
!>           (                  c(k)                    s(k) )
!>           (                         1                     )
!>           (                              ...              )
!>           (                                     1         )
!>           (                 -s(k)                    c(k) )
!> 
!> where R(k) appears in rows and columns k and z.  The rotations are
!> performed without ever forming P(k) explicitly.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          Specifies whether the plane rotation matrix P is applied to
!>          A on the left or the right.
!>          = 'L':  Left, compute A := P*A
!>          = 'R':  Right, compute A:= A*P**T
!> \endverbatim
!>
!> \param[in] PIVOT
!> \verbatim
!>          PIVOT is CHARACTER*1
!>          Specifies the plane for which P(k) is a plane rotation
!>          matrix.
!>          = 'V':  Variable pivot, the plane (k,k+1)
!>          = 'T':  Top pivot, the plane (1,k+1)
!>          = 'B':  Bottom pivot, the plane (k,z)
!> \endverbatim
!>
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Specifies whether P is a forward or backward sequence of
!>          plane rotations.
!>          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)
!>          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  If m <= 1, an immediate
!>          return is effected.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  If n <= 1, an
!>          immediate return is effected.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension
!>                  (M-1) if SIDE = 'L'
!>                  (N-1) if SIDE = 'R'
!>          The cosines c(k) of the plane rotations.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension
!>                  (M-1) if SIDE = 'L'
!>                  (N-1) if SIDE = 'R'
!>          The sines s(k) of the plane rotations.  The 2-by-2 plane
!>          rotation part of the matrix P(k), R(k), has the form
!>          R(k) = (  c(k)  s(k) )
!>                 ( -s(k)  c(k) ).
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The M-by-N matrix A.  On exit, A is overwritten by P*A if
!>          SIDE = 'R' or by A*P**T if SIDE = 'L'.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
    SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP, TEMP
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,&
           'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )&
           THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASR ', INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  P * A
!
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
10                   CONTINUE
                  END IF
20             CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!        Form A * P**T
!
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of DLASR
!
   END SUBROUTINE DLASR
!
!=
!
!> \brief \b DLASRT sorts numbers in increasing or decreasing order.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLASRT + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasrt.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasrt.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasrt.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASRT( ID, N, D, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          ID
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Sort the numbers in D in increasing order (if ID = 'I') or
!> in decreasing order (if ID = 'D' ).
!>
!> Use Quick Sort, reverting to Insertion sort on arrays of
!> size <= 20. Dimension of STACK limits N to about 2**32.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ID
!> \verbatim
!>          ID is CHARACTER*1
!>          = 'I': sort D in increasing order;
!>          = 'D': sort D in decreasing order.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The length of the array D.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          On entry, the array to be sorted.
!>          On exit, D has been sorted into increasing order
!>          (D(1) <= ... <= D(N) ) or into decreasing order
!>          (D(1) >= ... >= D(N) ), depending on ID.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
   SUBROUTINE DLASRT( ID, N, D, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
     CHARACTER          ID
     INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
     DOUBLE PRECISION   D( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
     INTEGER            SELECT
     PARAMETER          ( SELECT = 20 )
!     ..
!     .. Local Scalars ..
     INTEGER            DIR, ENDD, I, J, START, STKPNT
     DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
!     ..
!     .. Local Arrays ..
     INTEGER            STACK( 2, 32 )
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input paramters.
!
     INFO = 0
     DIR = -1
     IF( LSAME( ID, 'D' ) ) THEN
        DIR = 0
     ELSE IF( LSAME( ID, 'I' ) ) THEN
        DIR = 1
     END IF
     IF( DIR.EQ.-1 ) THEN
        INFO = -1
     ELSE IF( N.LT.0 ) THEN
        INFO = -2
     END IF
     IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLASRT', -INFO )
        RETURN
     END IF
!
!     Quick return if possible
!
     IF( N.LE.1 ) RETURN
!
     STKPNT = 1
     STACK( 1, 1 ) = 1
     STACK( 2, 1 ) = N
10   CONTINUE
     START = STACK( 1, STKPNT )
     ENDD = STACK( 2, STKPNT )
     STKPNT = STKPNT - 1
     IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
!
!        Do Insertion sort on D( START:ENDD )
!
        IF( DIR.EQ.0 ) THEN
!
!           Sort into decreasing order
!
           DO 30 I = START + 1, ENDD
              DO 20 J = I, START + 1, -1
                 IF( D( J ).GT.D( J-1 ) ) THEN
                    DMNMX = D( J )
                    D( J ) = D( J-1 )
                    D( J-1 ) = DMNMX
                 ELSE
                    GO TO 30
                 END IF
20           CONTINUE
30        CONTINUE
!
        ELSE
!
!           Sort into increasing order
!
           DO 50 I = START + 1, ENDD
              DO 40 J = I, START + 1, -1
                 IF( D( J ).LT.D( J-1 ) ) THEN
                    DMNMX = D( J )
                    D( J ) = D( J-1 )
                    D( J-1 ) = DMNMX
                 ELSE
                    GO TO 50
                 END IF
40            CONTINUE
50         CONTINUE
!
        END IF
!
     ELSE IF( ENDD-START.GT.SELECT ) THEN
!
!        Partition D( START:ENDD ) and stack parts, largest one first
!
!        Choose partition entry as median of 3
!
        D1 = D( START )
        D2 = D( ENDD )
        I = ( START+ENDD ) / 2
        D3 = D( I )
        IF( D1.LT.D2 ) THEN
           IF( D3.LT.D1 ) THEN
              DMNMX = D1
           ELSE IF( D3.LT.D2 ) THEN
              DMNMX = D3
           ELSE
              DMNMX = D2
           END IF
        ELSE
           IF( D3.LT.D2 ) THEN
              DMNMX = D2
           ELSE IF( D3.LT.D1 ) THEN
              DMNMX = D3
           ELSE
              DMNMX = D1
           END IF
        END IF
!
        IF( DIR.EQ.0 ) THEN
!
!           Sort into decreasing order
!
           I = START - 1
           J = ENDD + 1
60         CONTINUE
70         CONTINUE
           J = J - 1
           IF( D( J ).LT.DMNMX ) GO TO 70
80         CONTINUE
           I = I + 1
           IF( D( I ).GT.DMNMX ) GO TO 80
           IF( I.LT.J ) THEN
              TMP = D( I )
              D( I ) = D( J )
              D( J ) = TMP
              GO TO 60
           END IF
           IF( J-START.GT.ENDD-J-1 ) THEN
              STKPNT = STKPNT + 1
              STACK( 1, STKPNT ) = START
              STACK( 2, STKPNT ) = J
              STKPNT = STKPNT + 1
              STACK( 1, STKPNT ) = J + 1
              STACK( 2, STKPNT ) = ENDD
           ELSE
              STKPNT = STKPNT + 1
              STACK( 1, STKPNT ) = J + 1
              STACK( 2, STKPNT ) = ENDD
              STKPNT = STKPNT + 1
              STACK( 1, STKPNT ) = START
              STACK( 2, STKPNT ) = J
           END IF
        ELSE
!
!           Sort into increasing order
!
           I = START - 1
           J = ENDD + 1
90         CONTINUE
100        CONTINUE
           J = J - 1
           IF( D( J ).GT.DMNMX ) GO TO 100
110        CONTINUE
           I = I + 1
           IF( D( I ).LT.DMNMX ) GO TO 110
           IF( I.LT.J ) THEN
              TMP = D( I )
              D( I ) = D( J )
              D( J ) = TMP
              GO TO 90
           END IF
           IF( J-START.GT.ENDD-J-1 ) THEN
              STKPNT = STKPNT + 1
              STACK( 1, STKPNT ) = START
              STACK( 2, STKPNT ) = J
              STKPNT = STKPNT + 1
              STACK( 1, STKPNT ) = J + 1
              STACK( 2, STKPNT ) = ENDD
           ELSE
              STKPNT = STKPNT + 1
              STACK( 1, STKPNT ) = J + 1
              STACK( 2, STKPNT ) = ENDD
              STKPNT = STKPNT + 1
              STACK( 1, STKPNT ) = START
              STACK( 2, STKPNT ) = J
           END IF
        END IF
     END IF
     IF( STKPNT.GT.0 ) GO TO 10
     RETURN
!
!     End of DLASRT
!
  END SUBROUTINE DLASRT
!
!=
!
!> \brief \b DLAE2 computes the eigenvalues of a 2-by-2 symmetric matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAE2 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlae2.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlae2.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlae2.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   A, B, C, RT1, RT2
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
!>    [  A   B  ]
!>    [  B   C  ].
!> On return, RT1 is the eigenvalue of larger absolute value, and RT2
!> is the eigenvalue of smaller absolute value.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION
!>          The (1,1) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION
!>          The (1,2) and (2,1) elements of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!>          The (2,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[out] RT1
!> \verbatim
!>          RT1 is DOUBLE PRECISION
!>          The eigenvalue of larger absolute value.
!> \endverbatim
!>
!> \param[out] RT2
!> \verbatim
!>          RT2 is DOUBLE PRECISION
!>          The eigenvalue of smaller absolute value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  RT1 is accurate to a few ulps barring over/underflow.
!>
!>  RT2 may be inaccurate if there is massive cancellation in the
!>  determinant A*C-B*B; higher precision or correctly rounded or
!>  correctly truncated arithmetic would be needed to compute RT2
!>  accurately in all cases.
!>
!>  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!>  Underflow is harmless if the input data is 0 or exceeds
!>     underflow_threshold / macheps.
!> \endverbatim
!>
!  =====================================================================
  SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, RT1, RT2
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
!     Compute the eigenvalues
!
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
!
!        Includes case AB=ADF=0
!
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
!
!        Includes case RT1 = RT2 = 0
!
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
!
!     End of DLAE2
!
    END SUBROUTINE DLAE2
!
!=
!
!> \brief \b DLARTG generates a plane rotation with real cosine and real sine.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLARTG + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlartg.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlartg.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlartg.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARTG( F, G, CS, SN, R )
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   CS, F, G, R, SN
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARTG generate a plane rotation so that
!>
!>    [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
!>    [ -SN  CS  ]     [ G ]     [ 0 ]
!>
!> This is a slower, more accurate version of the BLAS1 routine DROTG,
!> with the following other differences:
!>    F and G are unchanged on return.
!>    If G=0, then CS=1 and SN=0.
!>    If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
!>       floating point operations (saves work in DBDSQR when
!>       there are zeros on the diagonal).
!>
!> If F exceeds G in magnitude, CS will be positive.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] F
!> \verbatim
!>          F is DOUBLE PRECISION
!>          The first component of vector to be rotated.
!> \endverbatim
!>
!> \param[in] G
!> \verbatim
!>          G is DOUBLE PRECISION
!>          The second component of vector to be rotated.
!> \endverbatim
!>
!> \param[out] CS
!> \verbatim
!>          CS is DOUBLE PRECISION
!>          The cosine of the rotation.
!> \endverbatim
!>
!> \param[out] SN
!> \verbatim
!>          SN is DOUBLE PRECISION
!>          The sine of the rotation.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is DOUBLE PRECISION
!>          The nonzero component of the rotated vector.
!>
!>  This version has a few statements commented out for thread safety
!>  (machine parameters are computed on each entry). 10 feb 03, SJH.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
    SUBROUTINE DLARTG( F, G, CS, SN, R )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
!     LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, INT, LOG, MAX, SQRT
!     ..
!     .. Save statement ..
!     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
!     ..
!     .. Data statements ..
!     DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
!     IF( FIRST ) THEN
      SAFMIN = DLAMCH( 'S' )
      EPS = DLAMCH( 'E' )
      SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /&
           LOG( DLAMCH( 'B' ) ) / TWO )
      SAFMX2 = ONE / SAFMN2
!        FIRST = .FALSE.
!     END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 ) GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 )  GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
!
!     End of DLARTG
!
   END SUBROUTINE DLARTG
!
!=
!
!> \brief \b DLACPY copies all or part of one two-dimensional array to another.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLACPY + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacpy.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacpy.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacpy.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDB, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLACPY copies all or part of a two-dimensional matrix A to another
!> matrix B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be copied to B.
!>          = 'U':      Upper triangular part
!>          = 'L':      Lower triangular part
!>          Otherwise:  All of the matrix A
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!>          or trapezoid is accessed; if UPLO = 'L', only the lower
!>          triangle or trapezoid is accessed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          On exit, B = A in the locations specified by UPLO.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
   SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
     CHARACTER          UPLO
     INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
     DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
     INTEGER            I, J
!     ..
!     .. External Functions ..
!     LOGICAL            LSAME
!     EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
     IF( LSAME( UPLO, 'U' ) ) THEN
        DO 20 J = 1, N
           DO 10 I = 1, MIN( J, M )
              B( I, J ) = A( I, J )
10         CONTINUE
20      CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
!
!     End of DLACPY
!
   END SUBROUTINE DLACPY
!
!=
!
!> \brief \b DLAED0 used by sstedc. Computes all eigenvalues and corresponding eigenvectors of an unreduced symmetric tridiagonal matrix using the divide and conquer method.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED0 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed0.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed0.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed0.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS,
!                          WORK, IWORK, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            ICOMPQ, INFO, LDQ, LDQS, N, QSIZ
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   D( * ), E( * ), Q( LDQ, * ), QSTORE( LDQS, * ),
!      $                   WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED0 computes all eigenvalues and corresponding eigenvectors of a
!> symmetric tridiagonal matrix using the divide and conquer method.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ICOMPQ
!> \verbatim
!>          ICOMPQ is INTEGER
!>          = 0:  Compute eigenvalues only.
!>          = 1:  Compute eigenvectors of original dense symmetric matrix
!>                also.  On entry, Q contains the orthogonal matrix used
!>                to reduce the original matrix to tridiagonal form.
!>          = 2:  Compute eigenvalues and eigenvectors of tridiagonal
!>                matrix.
!> \endverbatim
!>
!> \param[in] QSIZ
!> \verbatim
!>          QSIZ is INTEGER
!>         The dimension of the orthogonal matrix used to reduce
!>         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>         On entry, the main diagonal of the tridiagonal matrix.
!>         On exit, its eigenvalues.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>         The off-diagonal elements of the tridiagonal matrix.
!>         On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
!>         On entry, Q must contain an N-by-N orthogonal matrix.
!>         If ICOMPQ = 0    Q is not referenced.
!>         If ICOMPQ = 1    On entry, Q is a subset of the columns of the
!>                          orthogonal matrix used to reduce the full
!>                          matrix to tridiagonal form corresponding to
!>                          the subset of the full matrix which is being
!>                          decomposed at this time.
!>         If ICOMPQ = 2    On entry, Q will be the identity matrix.
!>                          On exit, Q contains the eigenvectors of the
!>                          tridiagonal matrix.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>         The leading dimension of the array Q.  If eigenvectors are
!>         desired, then  LDQ >= max(1,N).  In any case,  LDQ >= 1.
!> \endverbatim
!>
!> \param[out] QSTORE
!> \verbatim
!>          QSTORE is DOUBLE PRECISION array, dimension (LDQS, N)
!>         Referenced only when ICOMPQ = 1.  Used to store parts of
!>         the eigenvector matrix when the updating matrix multiplies
!>         take place.
!> \endverbatim
!>
!> \param[in] LDQS
!> \verbatim
!>          LDQS is INTEGER
!>         The leading dimension of the array QSTORE.  If ICOMPQ = 1,
!>         then  LDQS >= max(1,N).  In any case,  LDQS >= 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array,
!>         If ICOMPQ = 0 or 1, the dimension of WORK must be at least
!>                     1 + 3*N + 2*N*lg N + 3*N**2
!>                     ( lg( N ) = smallest integer k
!>                                 such that 2^k >= N )
!>         If ICOMPQ = 2, the dimension of WORK must be at least
!>                     4*N + N**2.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array,
!>         If ICOMPQ = 0 or 1, the dimension of IWORK must be at least
!>                        6 + 6*N + 5*N*lg N.
!>                        ( lg( N ) = smallest integer k
!>                                    such that 2^k >= N )
!>         If ICOMPQ = 2, the dimension of IWORK must be at least
!>                        3 + 5*N.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  The algorithm failed to compute an eigenvalue while
!>                working on the submatrix lying in rows and columns
!>                INFO/(N+1) through mod(INFO,N+1).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
   SUBROUTINE DLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS,&
        WORK, IWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
     INTEGER            ICOMPQ, INFO, LDQ, LDQS, N, QSIZ
!     ..
!     .. Array Arguments ..
     INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), Q( LDQ, * ), QSTORE( LDQS, * ),&
           WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            CURLVL, CURPRB, CURR, I, IGIVCL, IGIVNM,&
           IGIVPT, INDXQ, IPERM, IPRMPT, IQ, IQPTR, IWREM,&
           J, K, LGN, MATSIZ, MSD2, SMLSIZ, SMM1, SPM1,&
           SPM2, SUBMAT, SUBPBS, TLVLS
      DOUBLE PRECISION   TEMP
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DGEMM, DLACPY, DLAED1, DLAED7, DSTEQR,
!     $                   XERBLA
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, DBLE, INT, LOG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.2 ) THEN
         INFO = -1
      ELSE IF( ( ICOMPQ.EQ.1 ) .AND. ( QSIZ.LT.MAX( 0, N ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDQS.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED0', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
      SMLSIZ = ILAENV( 9, 'DLAED0', ' ', 0, 0, 0, 0 )
!
!     Determine the size and placement of the submatrices, and save in
!     the leading elements of IWORK.
!
      IWORK( 1 ) = N
      SUBPBS = 1
      TLVLS = 0
   10 CONTINUE
      IF( IWORK( SUBPBS ).GT.SMLSIZ ) THEN
         DO 20 J = SUBPBS, 1, -1
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2
            IWORK( 2*J-1 ) = IWORK( J ) / 2
   20    CONTINUE
         TLVLS = TLVLS + 1
         SUBPBS = 2*SUBPBS
         GO TO 10
      END IF
      DO 30 J = 2, SUBPBS
         IWORK( J ) = IWORK( J ) + IWORK( J-1 )
   30 CONTINUE
!
!     Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
!     using rank-1 modifications (cuts).
!
      SPM1 = SUBPBS - 1
      DO 40 I = 1, SPM1
         SUBMAT = IWORK( I ) + 1
         SMM1 = SUBMAT - 1
         D( SMM1 ) = D( SMM1 ) - ABS( E( SMM1 ) )
         D( SUBMAT ) = D( SUBMAT ) - ABS( E( SMM1 ) )
   40 CONTINUE
!
      INDXQ = 4*N + 3
      IF( ICOMPQ.NE.2 ) THEN
!
!        Set up workspaces for eigenvalues only/accumulate new vectors
!        routine
!
         TEMP = LOG( DBLE( N ) ) / LOG( TWO )
         LGN = INT( TEMP )
         IF( 2**LGN.LT.N ) LGN = LGN + 1
         IF( 2**LGN.LT.N ) LGN = LGN + 1
         IPRMPT = INDXQ + N + 1
         IPERM = IPRMPT + N*LGN
         IQPTR = IPERM + N*LGN
         IGIVPT = IQPTR + N + 2
         IGIVCL = IGIVPT + N*LGN
!
         IGIVNM = 1
         IQ = IGIVNM + 2*N*LGN
         IWREM = IQ + N**2 + 1
!
!        Initialize pointers
!
         DO 50 I = 0, SUBPBS
            IWORK( IPRMPT+I ) = 1
            IWORK( IGIVPT+I ) = 1
   50    CONTINUE
         IWORK( IQPTR ) = 1
      END IF
!
!     Solve each submatrix eigenproblem at the bottom of the divide and
!     conquer tree.
!
      CURR = 0
      DO 70 I = 0, SPM1
         IF( I.EQ.0 ) THEN
            SUBMAT = 1
            MATSIZ = IWORK( 1 )
         ELSE
            SUBMAT = IWORK( I ) + 1
            MATSIZ = IWORK( I+1 ) - IWORK( I )
         END IF
         IF( ICOMPQ.EQ.2 ) THEN
            CALL DSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ),&
                 Q( SUBMAT, SUBMAT ), LDQ, WORK, INFO )
            IF( INFO.NE.0 ) GO TO 130
         ELSE
            CALL DSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ),&
                 WORK( IQ-1+IWORK( IQPTR+CURR ) ), MATSIZ, WORK, INFO )
            IF( INFO.NE.0 ) GO TO 130
            IF( ICOMPQ.EQ.1 ) THEN
               CALL DGEMM( 'N', 'N', QSIZ, MATSIZ, MATSIZ, ONE,&
                    Q( 1, SUBMAT ), LDQ, WORK( IQ-1+IWORK( IQPTR+&
                    CURR ) ), MATSIZ, ZERO, QSTORE( 1, SUBMAT ),&
                    LDQS )
            END IF
            IWORK( IQPTR+CURR+1 ) = IWORK( IQPTR+CURR ) + MATSIZ**2
            CURR = CURR + 1
         END IF
         K = 1
         DO 60 J = SUBMAT, IWORK( I+1 )
            IWORK( INDXQ+J ) = K
            K = K + 1
   60    CONTINUE
   70 CONTINUE
!
!     Successively merge eigensystems of adjacent submatrices
!     into eigensystem for the corresponding larger matrix.
!
!     while ( SUBPBS > 1 )
!
      CURLVL = 1
   80 CONTINUE
      IF( SUBPBS.GT.1 ) THEN
         SPM2 = SUBPBS - 2
         DO 90 I = 0, SPM2, 2
            IF( I.EQ.0 ) THEN
               SUBMAT = 1
               MATSIZ = IWORK( 2 )
               MSD2 = IWORK( 1 )
               CURPRB = 0
            ELSE
               SUBMAT = IWORK( I ) + 1
               MATSIZ = IWORK( I+2 ) - IWORK( I )
               MSD2 = MATSIZ / 2
               CURPRB = CURPRB + 1
            END IF
!
!     Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
!     into an eigensystem of size MATSIZ.
!     DLAED1 is used only for the full eigensystem of a tridiagonal
!     matrix.
!     DLAED7 handles the cases in which eigenvalues only or eigenvalues
!     and eigenvectors of a full symmetric matrix (which was reduced to
!     tridiagonal form) are desired.
!
            IF( ICOMPQ.EQ.2 ) THEN
               CALL DLAED1( MATSIZ, D( SUBMAT ), Q( SUBMAT, SUBMAT ),&
                    LDQ, IWORK( INDXQ+SUBMAT ),&
                    E( SUBMAT+MSD2-1 ), MSD2, WORK,&
                    IWORK( SUBPBS+1 ), INFO )
            ELSE
               CALL DLAED7( ICOMPQ, MATSIZ, QSIZ, TLVLS, CURLVL, CURPRB,&
                    D( SUBMAT ), QSTORE( 1, SUBMAT ), LDQS,&
                    IWORK( INDXQ+SUBMAT ), E( SUBMAT+MSD2-1 ),&
                    MSD2, WORK( IQ ), IWORK( IQPTR ),&
                    IWORK( IPRMPT ), IWORK( IPERM ),&
                    IWORK( IGIVPT ), IWORK( IGIVCL ),&
                    WORK( IGIVNM ), WORK( IWREM ),&
                    IWORK( SUBPBS+1 ), INFO )
            END IF
            IF( INFO.NE.0 ) GO TO 130
            IWORK( I / 2+1 ) = IWORK( I+2 )
   90    CONTINUE
         SUBPBS = SUBPBS / 2
         CURLVL = CURLVL + 1
         GO TO 80
      END IF
!
!     end while
!
!     Re-merge the eigenvalues/vectors which were deflated at the final
!     merge step.
!
      IF( ICOMPQ.EQ.1 ) THEN
         DO 100 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
            CALL DCOPY( QSIZ, QSTORE( 1, J ), 1, Q( 1, I ), 1 )
  100    CONTINUE
         CALL DCOPY( N, WORK, 1, D, 1 )
      ELSE IF( ICOMPQ.EQ.2 ) THEN
         DO 110 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
            CALL DCOPY( N, Q( 1, J ), 1, WORK( N*I+1 ), 1 )
  110    CONTINUE
         CALL DCOPY( N, WORK, 1, D, 1 )
         CALL DLACPY( 'A', N, N, WORK( N+1 ), N, Q, LDQ )
      ELSE
         DO 120 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
  120    CONTINUE
         CALL DCOPY( N, WORK, 1, D, 1 )
      END IF
      GO TO 140
!
  130 CONTINUE
      INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1
!
  140 CONTINUE
      RETURN
!
!     End of DLAED0
!
   END SUBROUTINE DLAED0
!
!=
!
!> \brief \b ILADLC scans a matrix for its last non-zero column.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ILADLC + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlc.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladlc.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladlc.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILADLC( M, N, A, LDA )
! 
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILADLC scans A for its last non-zero column.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
   INTEGER FUNCTION ILADLC( M, N, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
     INTEGER            M, N, LDA
!     ..
!     .. Array Arguments ..
     DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
     DOUBLE PRECISION ZERO
     PARAMETER ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
     INTEGER I
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
     IF( N.EQ.0 ) THEN
        ILADLC = N
     ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
        ILADLC = N
     ELSE
!     Now scan each column from the end, returning with the first non-zero.
        DO ILADLC = N, 1, -1
           DO I = 1, M
              IF( A(I, ILADLC).NE.ZERO ) RETURN
           END DO
        END DO
      END IF
      RETURN
    END FUNCTION ILADLC
!
!=
!
!> \brief \b ILADLR scans a matrix for its last non-zero row.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ILADLR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladlr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladlr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILADLR( M, N, A, LDA )
! 
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILADLR scans A for its last non-zero row.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
    INTEGER FUNCTION ILADLR( M, N, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            M, N, LDA
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER I, J
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ILADLR = M
      ELSE IF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILADLR = M
      ELSE
!     Scan up each column tracking the last zero row seen.
         ILADLR = 0
         DO J = 1, N
            I=M
            DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
               I=I-1
            ENDDO
            ILADLR = MAX( ILADLR, I )
         END DO
      END IF
      RETURN
    END FUNCTION ILADLR
!
!=
!
!> \brief \b DLAED1 used by sstedc. Computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix. Used when the original matrix is tridiagonal.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED1 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed1.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed1.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed1.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED1( N, D, Q, LDQ, INDXQ, RHO, CUTPNT, WORK, IWORK,
!                          INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            CUTPNT, INFO, LDQ, N
!       DOUBLE PRECISION   RHO
!       ..
!       .. Array Arguments ..
!       INTEGER            INDXQ( * ), IWORK( * )
!       DOUBLE PRECISION   D( * ), Q( LDQ, * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED1 computes the updated eigensystem of a diagonal
!> matrix after modification by a rank-one symmetric matrix.  This
!> routine is used only for the eigenproblem which requires all
!> eigenvalues and eigenvectors of a tridiagonal matrix.  DLAED7 handles
!> the case in which eigenvalues only or eigenvalues and eigenvectors
!> of a full symmetric matrix (which was reduced to tridiagonal form)
!> are desired.
!>
!>   T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
!>
!>    where Z = Q**T*u, u is a vector of length N with ones in the
!>    CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
!>
!>    The eigenvectors of the original matrix are stored in Q, and the
!>    eigenvalues are in D.  The algorithm consists of three stages:
!>
!>       The first stage consists of deflating the size of the problem
!>       when there are multiple eigenvalues or if there is a zero in
!>       the Z vector.  For each such occurence the dimension of the
!>       secular equation problem is reduced by one.  This stage is
!>       performed by the routine DLAED2.
!>
!>       The second stage consists of calculating the updated
!>       eigenvalues. This is done by finding the roots of the secular
!>       equation via the routine DLAED4 (as called by DLAED3).
!>       This routine also calculates the eigenvectors of the current
!>       problem.
!>
!>       The final stage consists of computing the updated eigenvectors
!>       directly using the updated eigenvalues.  The eigenvectors for
!>       the current problem are multiplied with the eigenvectors from
!>       the overall problem.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>         On entry, the eigenvalues of the rank-1-perturbed matrix.
!>         On exit, the eigenvalues of the repaired matrix.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!>         On entry, the eigenvectors of the rank-1-perturbed matrix.
!>         On exit, the eigenvectors of the repaired tridiagonal matrix.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>         The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] INDXQ
!> \verbatim
!>          INDXQ is INTEGER array, dimension (N)
!>         On entry, the permutation which separately sorts the two
!>         subproblems in D into ascending order.
!>         On exit, the permutation which will reintegrate the
!>         subproblems back into sorted order,
!>         i.e. D( INDXQ( I = 1, N ) ) will be in ascending order.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>         The subdiagonal entry used to create the rank-1 modification.
!> \endverbatim
!>
!> \param[in] CUTPNT
!> \verbatim
!>          CUTPNT is INTEGER
!>         The location of the last eigenvalue in the leading sub-matrix.
!>         min(1,N) <= CUTPNT <= N/2.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (4*N + N**2)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (4*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, an eigenvalue did not converge
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA \n
!>  Modified by Francoise Tisseur, University of Tennessee
!>
!  =====================================================================
    SUBROUTINE DLAED1( N, D, Q, LDQ, INDXQ, RHO, CUTPNT, WORK, IWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            CUTPNT, INFO, LDQ, N
      DOUBLE PRECISION   RHO
!     ..
!     .. Array Arguments ..
      INTEGER            INDXQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), Q( LDQ, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            COLTYP, I, IDLMDA, INDX, INDXC, INDXP, IQ2, IS,&
           IW, IZ, K, N1, N2, ZPP1
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DLAED2, DLAED3, DLAMRG, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( MIN( 1, N / 2 ).GT.CUTPNT .OR. ( N / 2 ).LT.CUTPNT ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED1', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
!     The following values are integer pointers which indicate
!     the portion of the workspace
!     used by a particular array in DLAED2 and DLAED3.
!
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IQ2 = IW + N
!
      INDX = 1
      INDXC = INDX + N
      COLTYP = INDXC + N
      INDXP = COLTYP + N
!
!
!     Form the z-vector which consists of the last row of Q_1 and the
!     first row of Q_2.
!
      CALL DCOPY( CUTPNT, Q( CUTPNT, 1 ), LDQ, WORK( IZ ), 1 )
      ZPP1 = CUTPNT + 1
      CALL DCOPY( N-CUTPNT, Q( ZPP1, ZPP1 ), LDQ, WORK( IZ+CUTPNT ), 1 )
!
!     Deflate eigenvalues.
!
      CALL DLAED2( K, N, CUTPNT, D, Q, LDQ, INDXQ, RHO, WORK( IZ ),&
           WORK( IDLMDA ), WORK( IW ), WORK( IQ2 ),&
           IWORK( INDX ), IWORK( INDXC ), IWORK( INDXP ),&
           IWORK( COLTYP ), INFO )
!
      IF( INFO.NE.0 ) GO TO 20
!
!     Solve Secular Equation.
!
      IF( K.NE.0 ) THEN
         IS = ( IWORK( COLTYP )+IWORK( COLTYP+1 ) )*CUTPNT +&
              ( IWORK( COLTYP+1 )+IWORK( COLTYP+2 ) )*( N-CUTPNT ) + IQ2
         CALL DLAED3( K, N, CUTPNT, D, Q, LDQ, RHO, WORK( IDLMDA ),&
              WORK( IQ2 ), IWORK( INDXC ), IWORK( COLTYP ),&
              WORK( IW ), WORK( IS ), INFO )
         IF( INFO.NE.0 ) GO TO 20
!
!     Prepare the INDXQ sorting permutation.
!
         N1 = K
         N2 = N - K
         CALL DLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         DO 10 I = 1, N
            INDXQ( I ) = I
   10    CONTINUE
      END IF
!
   20 CONTINUE
      RETURN
!
!     End of DLAED1
!
   END SUBROUTINE DLAED1
!
!=
!
!> \brief \b DLAED7 used by sstedc. Computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix. Used when the original matrix is dense.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED7 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed7.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed7.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed7.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED7( ICOMPQ, N, QSIZ, TLVLS, CURLVL, CURPBM, D, Q,
!                          LDQ, INDXQ, RHO, CUTPNT, QSTORE, QPTR, PRMPTR,
!                          PERM, GIVPTR, GIVCOL, GIVNUM, WORK, IWORK,
!                          INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            CURLVL, CURPBM, CUTPNT, ICOMPQ, INFO, LDQ, N,
!      $                   QSIZ, TLVLS
!       DOUBLE PRECISION   RHO
!       ..
!       .. Array Arguments ..
!       INTEGER            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ),
!      $                   IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )
!       DOUBLE PRECISION   D( * ), GIVNUM( 2, * ), Q( LDQ, * ),
!      $                   QSTORE( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED7 computes the updated eigensystem of a diagonal
!> matrix after modification by a rank-one symmetric matrix. This
!> routine is used only for the eigenproblem which requires all
!> eigenvalues and optionally eigenvectors of a dense symmetric matrix
!> that has been reduced to tridiagonal form.  DLAED1 handles
!> the case in which all eigenvalues and eigenvectors of a symmetric
!> tridiagonal matrix are desired.
!>
!>   T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
!>
!>    where Z = Q**Tu, u is a vector of length N with ones in the
!>    CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
!>
!>    The eigenvectors of the original matrix are stored in Q, and the
!>    eigenvalues are in D.  The algorithm consists of three stages:
!>
!>       The first stage consists of deflating the size of the problem
!>       when there are multiple eigenvalues or if there is a zero in
!>       the Z vector.  For each such occurence the dimension of the
!>       secular equation problem is reduced by one.  This stage is
!>       performed by the routine DLAED8.
!>
!>       The second stage consists of calculating the updated
!>       eigenvalues. This is done by finding the roots of the secular
!>       equation via the routine DLAED4 (as called by DLAED9).
!>       This routine also calculates the eigenvectors of the current
!>       problem.
!>
!>       The final stage consists of computing the updated eigenvectors
!>       directly using the updated eigenvalues.  The eigenvectors for
!>       the current problem are multiplied with the eigenvectors from
!>       the overall problem.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ICOMPQ
!> \verbatim
!>          ICOMPQ is INTEGER
!>          = 0:  Compute eigenvalues only.
!>          = 1:  Compute eigenvectors of original dense symmetric matrix
!>                also.  On entry, Q contains the orthogonal matrix used
!>                to reduce the original matrix to tridiagonal form.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] QSIZ
!> \verbatim
!>          QSIZ is INTEGER
!>         The dimension of the orthogonal matrix used to reduce
!>         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
!> \endverbatim
!>
!> \param[in] TLVLS
!> \verbatim
!>          TLVLS is INTEGER
!>         The total number of merging levels in the overall divide and
!>         conquer tree.
!> \endverbatim
!>
!> \param[in] CURLVL
!> \verbatim
!>          CURLVL is INTEGER
!>         The current level in the overall merge routine,
!>         0 <= CURLVL <= TLVLS.
!> \endverbatim
!>
!> \param[in] CURPBM
!> \verbatim
!>          CURPBM is INTEGER
!>         The current problem in the current level in the overall
!>         merge routine (counting from upper left to lower right).
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>         On entry, the eigenvalues of the rank-1-perturbed matrix.
!>         On exit, the eigenvalues of the repaired matrix.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
!>         On entry, the eigenvectors of the rank-1-perturbed matrix.
!>         On exit, the eigenvectors of the repaired tridiagonal matrix.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>         The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[out] INDXQ
!> \verbatim
!>          INDXQ is INTEGER array, dimension (N)
!>         The permutation which will reintegrate the subproblem just
!>         solved back into sorted order, i.e., D( INDXQ( I = 1, N ) )
!>         will be in ascending order.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>         The subdiagonal element used to create the rank-1
!>         modification.
!> \endverbatim
!>
!> \param[in] CUTPNT
!> \verbatim
!>          CUTPNT is INTEGER
!>         Contains the location of the last eigenvalue in the leading
!>         sub-matrix.  min(1,N) <= CUTPNT <= N.
!> \endverbatim
!>
!> \param[in,out] QSTORE
!> \verbatim
!>          QSTORE is DOUBLE PRECISION array, dimension (N**2+1)
!>         Stores eigenvectors of submatrices encountered during
!>         divide and conquer, packed together. QPTR points to
!>         beginning of the submatrices.
!> \endverbatim
!>
!> \param[in,out] QPTR
!> \verbatim
!>          QPTR is INTEGER array, dimension (N+2)
!>         List of indices pointing to beginning of submatrices stored
!>         in QSTORE. The submatrices are numbered starting at the
!>         bottom left of the divide and conquer tree, from left to
!>         right and bottom to top.
!> \endverbatim
!>
!> \param[in] PRMPTR
!> \verbatim
!>          PRMPTR is INTEGER array, dimension (N lg N)
!>         Contains a list of pointers which indicate where in PERM a
!>         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
!>         indicates the size of the permutation and also the size of
!>         the full, non-deflated problem.
!> \endverbatim
!>
!> \param[in] PERM
!> \verbatim
!>          PERM is INTEGER array, dimension (N lg N)
!>         Contains the permutations (from deflation and sorting) to be
!>         applied to each eigenblock.
!> \endverbatim
!>
!> \param[in] GIVPTR
!> \verbatim
!>          GIVPTR is INTEGER array, dimension (N lg N)
!>         Contains a list of pointers which indicate where in GIVCOL a
!>         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
!>         indicates the number of Givens rotations.
!> \endverbatim
!>
!> \param[in] GIVCOL
!> \verbatim
!>          GIVCOL is INTEGER array, dimension (2, N lg N)
!>         Each pair of numbers indicates a pair of columns to take place
!>         in a Givens rotation.
!> \endverbatim
!>
!> \param[in] GIVNUM
!> \verbatim
!>          GIVNUM is DOUBLE PRECISION array, dimension (2, N lg N)
!>         Each number indicates the S value to be used in the
!>         corresponding Givens rotation.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (3*N+2*QSIZ*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (4*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, an eigenvalue did not converge
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2015
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
   SUBROUTINE DLAED7( ICOMPQ, N, QSIZ, TLVLS, CURLVL, CURPBM, D, Q,&
        LDQ, INDXQ, RHO, CUTPNT, QSTORE, QPTR, PRMPTR,&
        PERM, GIVPTR, GIVCOL, GIVNUM, WORK, IWORK,&
        INFO )
!
!  -- LAPACK computational routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
     INTEGER            CURLVL, CURPBM, CUTPNT, ICOMPQ, INFO, LDQ, N,&
          QSIZ, TLVLS
     DOUBLE PRECISION   RHO
!     ..
!     .. Array Arguments ..
     INTEGER            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ),&
          IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )
     DOUBLE PRECISION   D( * ), GIVNUM( 2, * ), Q( LDQ, * ),&
          QSTORE( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
     DOUBLE PRECISION   ONE, ZERO
     PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
     INTEGER            COLTYP, CURR, I, IDLMDA, INDX, INDXC, INDXP,&
          IQ2, IS, IW, IZ, K, LDQ2, N1, N2, PTR
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEMM, DLAED8, DLAED9, DLAEDA, DLAMRG, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
     INFO = 0
!
     IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
        INFO = -1
     ELSE IF( N.LT.0 ) THEN
        INFO = -2
     ELSE IF( ICOMPQ.EQ.1 .AND. QSIZ.LT.N ) THEN
        INFO = -3
     ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
        INFO = -9
     ELSE IF( MIN( 1, N ).GT.CUTPNT .OR. N.LT.CUTPNT ) THEN
        INFO = -12
     END IF
     IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLAED7', -INFO )
        RETURN
     END IF
!
!     Quick return if possible
!
     IF( N.EQ.0 ) RETURN
!
!     The following values are for bookkeeping purposes only.  They are
!     integer pointers which indicate the portion of the workspace
!     used by a particular array in DLAED8 and DLAED9.
!
     IF( ICOMPQ.EQ.1 ) THEN
        LDQ2 = QSIZ
     ELSE
        LDQ2 = N
     END IF
!
     IZ = 1
     IDLMDA = IZ + N
     IW = IDLMDA + N
     IQ2 = IW + N
     IS = IQ2 + N*LDQ2
!
     INDX = 1
     INDXC = INDX + N
     COLTYP = INDXC + N
     INDXP = COLTYP + N
!
!     Form the z-vector which consists of the last row of Q_1 and the
!     first row of Q_2.
!
     PTR = 1 + 2**TLVLS
     DO 10 I = 1, CURLVL - 1
        PTR = PTR + 2**( TLVLS-I )
10   CONTINUE
     CURR = PTR + CURPBM
     CALL DLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR,&
          GIVCOL, GIVNUM, QSTORE, QPTR, WORK( IZ ),&
          WORK( IZ+N ), INFO )
!
!     When solving the final problem, we no longer need the stored data,
!     so we will overwrite the data from this level onto the previously
!     used storage space.
!
     IF( CURLVL.EQ.TLVLS ) THEN
        QPTR( CURR ) = 1
        PRMPTR( CURR ) = 1
        GIVPTR( CURR ) = 1
     END IF
!
!     Sort and Deflate eigenvalues.
!
     CALL DLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO, CUTPNT,&
          WORK( IZ ), WORK( IDLMDA ), WORK( IQ2 ), LDQ2,&
          WORK( IW ), PERM( PRMPTR( CURR ) ), GIVPTR( CURR+1 ),&
          GIVCOL( 1, GIVPTR( CURR ) ),&
          GIVNUM( 1, GIVPTR( CURR ) ), IWORK( INDXP ),&
          IWORK( INDX ), INFO )
     PRMPTR( CURR+1 ) = PRMPTR( CURR ) + N
     GIVPTR( CURR+1 ) = GIVPTR( CURR+1 ) + GIVPTR( CURR )
!
!     Solve Secular Equation.
!
     IF( K.NE.0 ) THEN
        CALL DLAED9( K, 1, K, N, D, WORK( IS ), K, RHO, WORK( IDLMDA ),&
             WORK( IW ), QSTORE( QPTR( CURR ) ), K, INFO )
        IF( INFO.NE.0 ) GO TO 30
        IF( ICOMPQ.EQ.1 ) THEN
           CALL DGEMM( 'N', 'N', QSIZ, K, K, ONE, WORK( IQ2 ), LDQ2,&
                QSTORE( QPTR( CURR ) ), K, ZERO, Q, LDQ )
        END IF
        QPTR( CURR+1 ) = QPTR( CURR ) + K**2
!
!     Prepare the INDXQ sorting permutation.
!
        N1 = K
        N2 = N - K
        CALL DLAMRG( N1, N2, D, 1, -1, INDXQ )
     ELSE
        QPTR( CURR+1 ) = QPTR( CURR )
        DO 20 I = 1, N
           INDXQ( I ) = I
20      CONTINUE
     END IF
!
30   CONTINUE
     RETURN
!
!     End of DLAED7
!
  END SUBROUTINE DLAED7
!
!=
!
!> \brief \b DLAED8 used by sstedc. Merges eigenvalues and deflates secular equation. Used when the original matrix is dense.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED8 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed8.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed8.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed8.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO,
!                          CUTPNT, Z, DLAMDA, Q2, LDQ2, W, PERM, GIVPTR,
!                          GIVCOL, GIVNUM, INDXP, INDX, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            CUTPNT, GIVPTR, ICOMPQ, INFO, K, LDQ, LDQ2, N,
!      $                   QSIZ
!       DOUBLE PRECISION   RHO
!       ..
!       .. Array Arguments ..
!       INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ),
!      $                   INDXQ( * ), PERM( * )
!       DOUBLE PRECISION   D( * ), DLAMDA( * ), GIVNUM( 2, * ),
!      $                   Q( LDQ, * ), Q2( LDQ2, * ), W( * ), Z( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED8 merges the two sets of eigenvalues together into a single
!> sorted set.  Then it tries to deflate the size of the problem.
!> There are two ways in which deflation can occur:  when two or more
!> eigenvalues are close together or if there is a tiny element in the
!> Z vector.  For each such occurrence the order of the related secular
!> equation problem is reduced by one.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ICOMPQ
!> \verbatim
!>          ICOMPQ is INTEGER
!>          = 0:  Compute eigenvalues only.
!>          = 1:  Compute eigenvectors of original dense symmetric matrix
!>                also.  On entry, Q contains the orthogonal matrix used
!>                to reduce the original matrix to tridiagonal form.
!> \endverbatim
!>
!> \param[out] K
!> \verbatim
!>          K is INTEGER
!>         The number of non-deflated eigenvalues, and the order of the
!>         related secular equation.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] QSIZ
!> \verbatim
!>          QSIZ is INTEGER
!>         The dimension of the orthogonal matrix used to reduce
!>         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>         On entry, the eigenvalues of the two submatrices to be
!>         combined.  On exit, the trailing (N-K) updated eigenvalues
!>         (those which were deflated) sorted into increasing order.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!>         If ICOMPQ = 0, Q is not referenced.  Otherwise,
!>         on entry, Q contains the eigenvectors of the partially solved
!>         system which has been previously updated in matrix
!>         multiplies with other partially solved eigensystems.
!>         On exit, Q contains the trailing (N-K) updated eigenvectors
!>         (those which were deflated) in its last N-K columns.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>         The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in] INDXQ
!> \verbatim
!>          INDXQ is INTEGER array, dimension (N)
!>         The permutation which separately sorts the two sub-problems
!>         in D into ascending order.  Note that elements in the second
!>         half of this permutation must first have CUTPNT added to
!>         their values in order to be accurate.
!> \endverbatim
!>
!> \param[in,out] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>         On entry, the off-diagonal element associated with the rank-1
!>         cut which originally split the two submatrices which are now
!>         being recombined.
!>         On exit, RHO has been modified to the value required by
!>         DLAED3.
!> \endverbatim
!>
!> \param[in] CUTPNT
!> \verbatim
!>          CUTPNT is INTEGER
!>         The location of the last eigenvalue in the leading
!>         sub-matrix.  min(1,N) <= CUTPNT <= N.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (N)
!>         On entry, Z contains the updating vector (the last row of
!>         the first sub-eigenvector matrix and the first row of the
!>         second sub-eigenvector matrix).
!>         On exit, the contents of Z are destroyed by the updating
!>         process.
!> \endverbatim
!>
!> \param[out] DLAMDA
!> \verbatim
!>          DLAMDA is DOUBLE PRECISION array, dimension (N)
!>         A copy of the first K eigenvalues which will be used by
!>         DLAED3 to form the secular equation.
!> \endverbatim
!>
!> \param[out] Q2
!> \verbatim
!>          Q2 is DOUBLE PRECISION array, dimension (LDQ2,N)
!>         If ICOMPQ = 0, Q2 is not referenced.  Otherwise,
!>         a copy of the first K eigenvectors which will be used by
!>         DLAED7 in a matrix multiply (DGEMM) to update the new
!>         eigenvectors.
!> \endverbatim
!>
!> \param[in] LDQ2
!> \verbatim
!>          LDQ2 is INTEGER
!>         The leading dimension of the array Q2.  LDQ2 >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>         The first k values of the final deflation-altered z-vector and
!>         will be passed to DLAED3.
!> \endverbatim
!>
!> \param[out] PERM
!> \verbatim
!>          PERM is INTEGER array, dimension (N)
!>         The permutations (from deflation and sorting) to be applied
!>         to each eigenblock.
!> \endverbatim
!>
!> \param[out] GIVPTR
!> \verbatim
!>          GIVPTR is INTEGER
!>         The number of Givens rotations which took place in this
!>         subproblem.
!> \endverbatim
!>
!> \param[out] GIVCOL
!> \verbatim
!>          GIVCOL is INTEGER array, dimension (2, N)
!>         Each pair of numbers indicates a pair of columns to take place
!>         in a Givens rotation.
!> \endverbatim
!>
!> \param[out] GIVNUM
!> \verbatim
!>          GIVNUM is DOUBLE PRECISION array, dimension (2, N)
!>         Each number indicates the S value to be used in the
!>         corresponding Givens rotation.
!> \endverbatim
!>
!> \param[out] INDXP
!> \verbatim
!>          INDXP is INTEGER array, dimension (N)
!>         The permutation used to place deflated values of D at the end
!>         of the array.  INDXP(1:K) points to the nondeflated D-values
!>         and INDXP(K+1:N) points to the deflated eigenvalues.
!> \endverbatim
!>
!> \param[out] INDX
!> \verbatim
!>          INDX is INTEGER array, dimension (N)
!>         The permutation used to sort the contents of D into ascending
!>         order.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
  SUBROUTINE DLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO,&
       CUTPNT, Z, DLAMDA, Q2, LDQ2, W, PERM, GIVPTR,&
       GIVCOL, GIVNUM, INDXP, INDX, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
    INTEGER            CUTPNT, GIVPTR, ICOMPQ, INFO, K, LDQ, LDQ2, N,&
         QSIZ
    DOUBLE PRECISION   RHO
!     ..
!     .. Array Arguments ..
    INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ),&
         INDXQ( * ), PERM( * )
    DOUBLE PRECISION   D( * ), DLAMDA( * ), GIVNUM( 2, * ),&
         Q( LDQ, * ), Q2( LDQ2, * ), W( * ), Z( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
    DOUBLE PRECISION   MONE, ZERO, ONE, TWO, EIGHT
    PARAMETER          ( MONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0,&
         TWO = 2.0D0, EIGHT = 8.0D0 )
!     ..
!     .. Local Scalars ..
!
    INTEGER            I, IMAX, J, JLAM, JMAX, JP, K2, N1, N1P1, N2
    DOUBLE PRECISION   C, EPS, S, T, TAU, TOL
!     ..
!     .. External Functions ..
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DLAMCH, DLAPY2
!      EXTERNAL           IDAMAX, DLAMCH, DLAPY2
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DLACPY, DLAMRG, DROT, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
    INFO = 0
!
    IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
       INFO = -1
    ELSE IF( N.LT.0 ) THEN
       INFO = -3
    ELSE IF( ICOMPQ.EQ.1 .AND. QSIZ.LT.N ) THEN
       INFO = -4
    ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
       INFO = -7
    ELSE IF( CUTPNT.LT.MIN( 1, N ) .OR. CUTPNT.GT.N ) THEN
       INFO = -10
    ELSE IF( LDQ2.LT.MAX( 1, N ) ) THEN
       INFO = -14
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DLAED8', -INFO )
       RETURN
    END IF
!
!     Need to initialize GIVPTR to O here in case of quick exit
!     to prevent an unspecified code behavior (usually sigfault) 
!     when IWORK array on entry to *stedc is not zeroed 
!     (or at least some IWORK entries which used in *laed7 for GIVPTR).
!
    GIVPTR = 0
!
!     Quick return if possible
!
    IF( N.EQ.0 ) RETURN
!
    N1 = CUTPNT
    N2 = N - N1
    N1P1 = N1 + 1
!
    IF( RHO.LT.ZERO ) THEN
       CALL DSCAL( N2, MONE, Z( N1P1 ), 1 )
    END IF
!
!     Normalize z so that norm(z) = 1
!
    T = ONE / SQRT( TWO )
    DO 10 J = 1, N
       INDX( J ) = J
10  CONTINUE
    CALL DSCAL( N, T, Z, 1 )
    RHO = ABS( TWO*RHO )
!
!     Sort the eigenvalues into increasing order
!
    DO 20 I = CUTPNT + 1, N
       INDXQ( I ) = INDXQ( I ) + CUTPNT
20  CONTINUE
    DO 30 I = 1, N
       DLAMDA( I ) = D( INDXQ( I ) )
       W( I ) = Z( INDXQ( I ) )
30  CONTINUE
    I = 1
    J = CUTPNT + 1
    CALL DLAMRG( N1, N2, DLAMDA, 1, 1, INDX )
    DO 40 I = 1, N
       D( I ) = DLAMDA( INDX( I ) )
       Z( I ) = W( INDX( I ) )
40  CONTINUE
!
!     Calculate the allowable deflation tolerence
!
    IMAX = IDAMAX( N, Z, 1 )
    JMAX = IDAMAX( N, D, 1 )
    EPS = DLAMCH( 'Epsilon' )
    TOL = EIGHT*EPS*ABS( D( JMAX ) )
!
!     If the rank-1 modifier is small enough, no more needs to be done
!     except to reorganize Q so that its columns correspond with the
!     elements in D.
!
    IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
       K = 0
       IF( ICOMPQ.EQ.0 ) THEN
          DO 50 J = 1, N
             PERM( J ) = INDXQ( INDX( J ) )
   50     CONTINUE
       ELSE
          DO 60 J = 1, N
             PERM( J ) = INDXQ( INDX( J ) )
             CALL DCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
   60     CONTINUE
          CALL DLACPY( 'A', QSIZ, N, Q2( 1, 1 ), LDQ2, Q( 1, 1 ), LDQ )
       END IF
       RETURN
    END IF
!
!     If there are multiple eigenvalues then the problem deflates.  Here
!     the number of equal eigenvalues are found.  As each equal
!     eigenvalue is found, an elementary reflector is computed to rotate
!     the corresponding eigensubspace so that the corresponding
!     components of Z are zero in this new basis.
!
    K = 0
    K2 = N + 1
    DO 70 J = 1, N
       IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
!
!           Deflate due to small z component.
!
          K2 = K2 - 1
          INDXP( K2 ) = J
          IF( J.EQ.N ) GO TO 110
       ELSE
          JLAM = J
          GO TO 80
       END IF
70  CONTINUE
80  CONTINUE
    J = J + 1
    IF( J.GT.N ) GO TO 100
    IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
!
!        Deflate due to small z component.
!
       K2 = K2 - 1
       INDXP( K2 ) = J
    ELSE
!
!        Check if eigenvalues are close enough to allow deflation.
!
       S = Z( JLAM )
       C = Z( J )
!
!        Find sqrt(a**2+b**2) without overflow or
!        destructive underflow.
!
       TAU = DLAPY2( C, S )
       T = D( J ) - D( JLAM )
       C = C / TAU
       S = -S / TAU
       IF( ABS( T*C*S ).LE.TOL ) THEN
!
!           Deflation is possible.
!
          Z( J ) = TAU
          Z( JLAM ) = ZERO
!
!           Record the appropriate Givens rotation
!
          GIVPTR = GIVPTR + 1
          GIVCOL( 1, GIVPTR ) = INDXQ( INDX( JLAM ) )
          GIVCOL( 2, GIVPTR ) = INDXQ( INDX( J ) )
          GIVNUM( 1, GIVPTR ) = C
          GIVNUM( 2, GIVPTR ) = S
          IF( ICOMPQ.EQ.1 ) THEN
             CALL DROT( QSIZ, Q( 1, INDXQ( INDX( JLAM ) ) ), 1,&
                  Q( 1, INDXQ( INDX( J ) ) ), 1, C, S )
          END IF
          T = D( JLAM )*C*C + D( J )*S*S
          D( J ) = D( JLAM )*S*S + D( J )*C*C
          D( JLAM ) = T
          K2 = K2 - 1
          I = 1
90        CONTINUE
          IF( K2+I.LE.N ) THEN
             IF( D( JLAM ).LT.D( INDXP( K2+I ) ) ) THEN
                INDXP( K2+I-1 ) = INDXP( K2+I )
                INDXP( K2+I ) = JLAM
                I = I + 1
                GO TO 90
             ELSE
                INDXP( K2+I-1 ) = JLAM
             END IF
          ELSE
             INDXP( K2+I-1 ) = JLAM
          END IF
          JLAM = J
       ELSE
          K = K + 1
          W( K ) = Z( JLAM )
          DLAMDA( K ) = D( JLAM )
          INDXP( K ) = JLAM
          JLAM = J
       END IF
    END IF
    GO TO 80
100 CONTINUE
!
!     Record the last eigenvalue.
!
    K = K + 1
    W( K ) = Z( JLAM )
    DLAMDA( K ) = D( JLAM )
    INDXP( K ) = JLAM
!
110 CONTINUE
!
!     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
!     and Q2 respectively.  The eigenvalues/vectors which were not
!     deflated go into the first K slots of DLAMDA and Q2 respectively,
!     while those which were deflated go into the last N - K slots.
!
    IF( ICOMPQ.EQ.0 ) THEN
       DO 120 J = 1, N
          JP = INDXP( J )
          DLAMDA( J ) = D( JP )
          PERM( J ) = INDXQ( INDX( JP ) )
120    CONTINUE
    ELSE
       DO 130 J = 1, N
          JP = INDXP( J )
          DLAMDA( J ) = D( JP )
          PERM( J ) = INDXQ( INDX( JP ) )
          CALL DCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
130    CONTINUE
    END IF
!
!     The deflated eigenvalues and their corresponding vectors go back
!     into the last N - K slots of D and Q respectively.
!
    IF( K.LT.N ) THEN
       IF( ICOMPQ.EQ.0 ) THEN
          CALL DCOPY( N-K, DLAMDA( K+1 ), 1, D( K+1 ), 1 )
       ELSE
          CALL DCOPY( N-K, DLAMDA( K+1 ), 1, D( K+1 ), 1 )
          CALL DLACPY( 'A', QSIZ, N-K, Q2( 1, K+1 ), LDQ2,&
               Q( 1, K+1 ), LDQ )
       END IF
    END IF
!
    RETURN
!
!     End of DLAED8
!
 END SUBROUTINE DLAED8
!
!=
!
!> \brief \b DLAED9 used by sstedc. Finds the roots of the secular equation and updates the eigenvectors. Used when the original matrix is dense.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED9 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed9.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed9.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed9.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMDA, W,
!                          S, LDS, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, KSTART, KSTOP, LDQ, LDS, N
!       DOUBLE PRECISION   RHO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), S( LDS, * ),
!      $                   W( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED9 finds the roots of the secular equation, as defined by the
!> values in D, Z, and RHO, between KSTART and KSTOP.  It makes the
!> appropriate calls to DLAED4 and then stores the new matrix of
!> eigenvectors for use in calculating the next level of Z vectors.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of terms in the rational function to be solved by
!>          DLAED4.  K >= 0.
!> \endverbatim
!>
!> \param[in] KSTART
!> \verbatim
!>          KSTART is INTEGER
!> \endverbatim
!>
!> \param[in] KSTOP
!> \verbatim
!>          KSTOP is INTEGER
!>          The updated eigenvalues Lambda(I), KSTART <= I <= KSTOP
!>          are to be computed.  1 <= KSTART <= KSTOP <= K.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns in the Q matrix.
!>          N >= K (delation may result in N > K).
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          D(I) contains the updated eigenvalues
!>          for KSTART <= I <= KSTOP.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max( 1, N ).
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>          The value of the parameter in the rank one update equation.
!>          RHO >= 0 required.
!> \endverbatim
!>
!> \param[in] DLAMDA
!> \verbatim
!>          DLAMDA is DOUBLE PRECISION array, dimension (K)
!>          The first K elements of this array contain the old roots
!>          of the deflated updating problem.  These are the poles
!>          of the secular equation.
!> \endverbatim
!>
!> \param[in] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (K)
!>          The first K elements of this array contain the components
!>          of the deflation-adjusted updating vector.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (LDS, K)
!>          Will contain the eigenvectors of the repaired matrix which
!>          will be stored for subsequent Z vector calculation and
!>          multiplied by the previously accumulated eigenvectors
!>          to update the system.
!> \endverbatim
!>
!> \param[in] LDS
!> \verbatim
!>          LDS is INTEGER
!>          The leading dimension of S.  LDS >= max( 1, K ).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, an eigenvalue did not converge
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
 SUBROUTINE DLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMDA, W,&
      S, LDS, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
   INTEGER            INFO, K, KSTART, KSTOP, LDQ, LDS, N
   DOUBLE PRECISION   RHO
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), S( LDS, * ),&
        W( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   INTEGER            I, J
   DOUBLE PRECISION   TEMP
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMC3, DNRM2
!      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DLAED4, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
!
   IF( K.LT.0 ) THEN
      INFO = -1
   ELSE IF( KSTART.LT.1 .OR. KSTART.GT.MAX( 1, K ) ) THEN
      INFO = -2
   ELSE IF( MAX( 1, KSTOP ).LT.KSTART .OR. KSTOP.GT.MAX( 1, K ) ) THEN
      INFO = -3
   ELSE IF( N.LT.K ) THEN
      INFO = -4
   ELSE IF( LDQ.LT.MAX( 1, K ) ) THEN
      INFO = -7
   ELSE IF( LDS.LT.MAX( 1, K ) ) THEN
      INFO = -12
   END IF
   IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'DLAED9', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( K.EQ.0 ) RETURN
!
!     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DLAMDA(I) if it is 1; this makes the subsequent
!     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DLAMDA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DLAMDA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
!     this code.
!
   DO 10 I = 1, N
      DLAMDA( I ) = DLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
10 CONTINUE
!
   DO 20 J = KSTART, KSTOP
      CALL DLAED4( K, J, DLAMDA, W, Q( 1, J ), RHO, D( J ), INFO )
!
!        If the zero finder fails, the computation is terminated.
!
      IF( INFO.NE.0 ) GO TO 120
20 CONTINUE
!
   IF( K.EQ.1 .OR. K.EQ.2 ) THEN
      DO 40 I = 1, K
         DO 30 J = 1, K
            S( J, I ) = Q( J, I )
30       CONTINUE
40    CONTINUE
      GO TO 120
   END IF
!
!     Compute updated W.
!
   CALL DCOPY( K, W, 1, S, 1 )
!
!     Initialize W(I) = Q(I,I)
!
   CALL DCOPY( K, Q, LDQ+1, W, 1 )
   DO 70 J = 1, K
      DO 50 I = 1, J - 1
         W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
50    CONTINUE
      DO 60 I = J + 1, K
         W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
60    CONTINUE
70 CONTINUE
   DO 80 I = 1, K
      W( I ) = SIGN( SQRT( -W( I ) ), S( I, 1 ) )
80 CONTINUE
!
!     Compute eigenvectors of the modified rank-1 modification.
!
   DO 110 J = 1, K
      DO 90 I = 1, K
         Q( I, J ) = W( I ) / Q( I, J )
90    CONTINUE
      TEMP = DNRM2( K, Q( 1, J ), 1 )
      DO 100 I = 1, K
         S( I, J ) = Q( I, J ) / TEMP
100   CONTINUE
110 CONTINUE
!
120 CONTINUE
   RETURN
!
!     End of DLAED9
!
 END SUBROUTINE DLAED9
!
!=
!
!> \brief \b DLAED2 used by sstedc. Merges eigenvalues and deflates secular equation. Used when the original matrix is tridiagonal.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED2 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed2.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed2.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed2.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED2( K, N, N1, D, Q, LDQ, INDXQ, RHO, Z, DLAMDA, W,
!                          Q2, INDX, INDXC, INDXP, COLTYP, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDQ, N, N1
!       DOUBLE PRECISION   RHO
!       ..
!       .. Array Arguments ..
!       INTEGER            COLTYP( * ), INDX( * ), INDXC( * ), INDXP( * ),
!      $                   INDXQ( * )
!       DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),
!      $                   W( * ), Z( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED2 merges the two sets of eigenvalues together into a single
!> sorted set.  Then it tries to deflate the size of the problem.
!> There are two ways in which deflation can occur:  when two or more
!> eigenvalues are close together or if there is a tiny entry in the
!> Z vector.  For each such occurrence the order of the related secular
!> equation problem is reduced by one.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] K
!> \verbatim
!>          K is INTEGER
!>         The number of non-deflated eigenvalues, and the order of the
!>         related secular equation. 0 <= K <=N.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>         The location of the last eigenvalue in the leading sub-matrix.
!>         min(1,N) <= N1 <= N/2.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>         On entry, D contains the eigenvalues of the two submatrices to
!>         be combined.
!>         On exit, D contains the trailing (N-K) updated eigenvalues
!>         (those which were deflated) sorted into increasing order.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
!>         On entry, Q contains the eigenvectors of two submatrices in
!>         the two square blocks with corners at (1,1), (N1,N1)
!>         and (N1+1, N1+1), (N,N).
!>         On exit, Q contains the trailing (N-K) updated eigenvectors
!>         (those which were deflated) in its last N-K columns.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>         The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] INDXQ
!> \verbatim
!>          INDXQ is INTEGER array, dimension (N)
!>         The permutation which separately sorts the two sub-problems
!>         in D into ascending order.  Note that elements in the second
!>         half of this permutation must first have N1 added to their
!>         values. Destroyed on exit.
!> \endverbatim
!>
!> \param[in,out] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>         On entry, the off-diagonal element associated with the rank-1
!>         cut which originally split the two submatrices which are now
!>         being recombined.
!>         On exit, RHO has been modified to the value required by
!>         DLAED3.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (N)
!>         On entry, Z contains the updating vector (the last
!>         row of the first sub-eigenvector matrix and the first row of
!>         the second sub-eigenvector matrix).
!>         On exit, the contents of Z have been destroyed by the updating
!>         process.
!> \endverbatim
!>
!> \param[out] DLAMDA
!> \verbatim
!>          DLAMDA is DOUBLE PRECISION array, dimension (N)
!>         A copy of the first K eigenvalues which will be used by
!>         DLAED3 to form the secular equation.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>         The first k values of the final deflation-altered z-vector
!>         which will be passed to DLAED3.
!> \endverbatim
!>
!> \param[out] Q2
!> \verbatim
!>          Q2 is DOUBLE PRECISION array, dimension (N1**2+(N-N1)**2)
!>         A copy of the first K eigenvectors which will be used by
!>         DLAED3 in a matrix multiply (DGEMM) to solve for the new
!>         eigenvectors.
!> \endverbatim
!>
!> \param[out] INDX
!> \verbatim
!>          INDX is INTEGER array, dimension (N)
!>         The permutation used to sort the contents of DLAMDA into
!>         ascending order.
!> \endverbatim
!>
!> \param[out] INDXC
!> \verbatim
!>          INDXC is INTEGER array, dimension (N)
!>         The permutation used to arrange the columns of the deflated
!>         Q matrix into three groups:  the first group contains non-zero
!>         elements only at and above N1, the second contains
!>         non-zero elements only below N1, and the third is dense.
!> \endverbatim
!>
!> \param[out] INDXP
!> \verbatim
!>          INDXP is INTEGER array, dimension (N)
!>         The permutation used to place deflated values of D at the end
!>         of the array.  INDXP(1:K) points to the nondeflated D-values
!>         and INDXP(K+1:N) points to the deflated eigenvalues.
!> \endverbatim
!>
!> \param[out] COLTYP
!> \verbatim
!>          COLTYP is INTEGER array, dimension (N)
!>         During execution, a label which will indicate which of the
!>         following types a column in the Q2 matrix is:
!>         1 : non-zero in the upper half only;
!>         2 : dense;
!>         3 : non-zero in the lower half only;
!>         4 : deflated.
!>         On exit, COLTYP(i) is the number of columns of type i,
!>         for i=1 to 4 only.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA \n
!>  Modified by Francoise Tisseur, University of Tennessee
!>
!  =====================================================================
 SUBROUTINE DLAED2( K, N, N1, D, Q, LDQ, INDXQ, RHO, Z, DLAMDA, W,&
      Q2, INDX, INDXC, INDXP, COLTYP, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
   INTEGER            INFO, K, LDQ, N, N1
   DOUBLE PRECISION   RHO
!     ..
!     .. Array Arguments ..
   INTEGER            COLTYP( * ), INDX( * ), INDXC( * ), INDXP( * ),&
        INDXQ( * )
   DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),&
        W( * ), Z( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   MONE, ZERO, ONE, TWO, EIGHT
   PARAMETER          ( MONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0,&
        TWO = 2.0D0, EIGHT = 8.0D0 )
!     ..
!     .. Local Arrays ..
   INTEGER            CTOT( 4 ), PSM( 4 )
!     ..
!     .. Local Scalars ..
   INTEGER            CT, I, IMAX, IQ1, IQ2, J, JMAX, JS, K2, N1P1,&
        N2, NJ, PJ
   DOUBLE PRECISION   C, EPS, S, T, TAU, TOL
!     ..
!     .. External Functions ..
!   INTEGER            IDAMAX
!   DOUBLE PRECISION   DLAMCH, DLAPY2
!   EXTERNAL           IDAMAX, DLAMCH, DLAPY2
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DLACPY, DLAMRG, DROT, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
!
   IF( N.LT.0 ) THEN
      INFO = -2
   ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
      INFO = -6
   ELSE IF( MIN( 1, ( N / 2 ) ).GT.N1 .OR. ( N / 2 ).LT.N1 ) THEN
      INFO = -3
   END IF
   IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'DLAED2', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( N.EQ.0 ) RETURN
!
   N2 = N - N1
   N1P1 = N1 + 1
!
   IF( RHO.LT.ZERO ) THEN
      CALL DSCAL( N2, MONE, Z( N1P1 ), 1 )
   END IF
!
!     Normalize z so that norm(z) = 1.  Since z is the concatenation of
!     two normalized vectors, norm2(z) = sqrt(2).
!
   T = ONE / SQRT( TWO )
   CALL DSCAL( N, T, Z, 1 )
!
!     RHO = ABS( norm(z)**2 * RHO )
!
   RHO = ABS( TWO*RHO )
!
!     Sort the eigenvalues into increasing order
!
   DO 10 I = N1P1, N
      INDXQ( I ) = INDXQ( I ) + N1
10 CONTINUE
!
!     re-integrate the deflated parts from the last pass
!
   DO 20 I = 1, N
      DLAMDA( I ) = D( INDXQ( I ) )
20 CONTINUE
   CALL DLAMRG( N1, N2, DLAMDA, 1, 1, INDXC )
   DO 30 I = 1, N
      INDX( I ) = INDXQ( INDXC( I ) )
30 CONTINUE
!
!     Calculate the allowable deflation tolerance
!
   IMAX = IDAMAX( N, Z, 1 )
   JMAX = IDAMAX( N, D, 1 )
   EPS = DLAMCH( 'Epsilon' )
   TOL = EIGHT*EPS*MAX( ABS( D( JMAX ) ), ABS( Z( IMAX ) ) )
!
!     If the rank-1 modifier is small enough, no more needs to be done
!     except to reorganize Q so that its columns correspond with the
!     elements in D.
!
   IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
      K = 0
      IQ2 = 1
      DO 40 J = 1, N
         I = INDX( J )
         CALL DCOPY( N, Q( 1, I ), 1, Q2( IQ2 ), 1 )
         DLAMDA( J ) = D( I )
         IQ2 = IQ2 + N
40    CONTINUE
      CALL DLACPY( 'A', N, N, Q2, N, Q, LDQ )
      CALL DCOPY( N, DLAMDA, 1, D, 1 )
      GO TO 190
   END IF
!
!     If there are multiple eigenvalues then the problem deflates.  Here
!     the number of equal eigenvalues are found.  As each equal
!     eigenvalue is found, an elementary reflector is computed to rotate
!     the corresponding eigensubspace so that the corresponding
!     components of Z are zero in this new basis.
!
   DO 50 I = 1, N1
      COLTYP( I ) = 1
50 CONTINUE
   DO 60 I = N1P1, N
      COLTYP( I ) = 3
60 CONTINUE
!
!
   K = 0
   K2 = N + 1
   DO 70 J = 1, N
      NJ = INDX( J )
      IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
!
!           Deflate due to small z component.
!
         K2 = K2 - 1
         COLTYP( NJ ) = 4
         INDXP( K2 ) = NJ
         IF( J.EQ.N ) GO TO 100
      ELSE
         PJ = NJ
         GO TO 80
      END IF
70 CONTINUE
80 CONTINUE
   J = J + 1
   NJ = INDX( J )
   IF( J.GT.N ) GO TO 100
   IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
!
!        Deflate due to small z component.
!
      K2 = K2 - 1
      COLTYP( NJ ) = 4
      INDXP( K2 ) = NJ
   ELSE
!
!        Check if eigenvalues are close enough to allow deflation.
!
      S = Z( PJ )
      C = Z( NJ )
!
!        Find sqrt(a**2+b**2) without overflow or
!        destructive underflow.
!
      TAU = DLAPY2( C, S )
      T = D( NJ ) - D( PJ )
      C = C / TAU
      S = -S / TAU
      IF( ABS( T*C*S ).LE.TOL ) THEN
!
!           Deflation is possible.
!
         Z( NJ ) = TAU
         Z( PJ ) = ZERO
         IF( COLTYP( NJ ).NE.COLTYP( PJ ) ) COLTYP( NJ ) = 2
         COLTYP( PJ ) = 4
         CALL DROT( N, Q( 1, PJ ), 1, Q( 1, NJ ), 1, C, S )
         T = D( PJ )*C**2 + D( NJ )*S**2
         D( NJ ) = D( PJ )*S**2 + D( NJ )*C**2
         D( PJ ) = T
         K2 = K2 - 1
         I = 1
90       CONTINUE
         IF( K2+I.LE.N ) THEN
            IF( D( PJ ).LT.D( INDXP( K2+I ) ) ) THEN
               INDXP( K2+I-1 ) = INDXP( K2+I )
               INDXP( K2+I ) = PJ
               I = I + 1
               GO TO 90
            ELSE
               INDXP( K2+I-1 ) = PJ
            END IF
         ELSE
            INDXP( K2+I-1 ) = PJ
         END IF
         PJ = NJ
      ELSE
         K = K + 1
         DLAMDA( K ) = D( PJ )
         W( K ) = Z( PJ )
         INDXP( K ) = PJ
         PJ = NJ
      END IF
   END IF
   GO TO 80
100 CONTINUE
!
!     Record the last eigenvalue.
!
   K = K + 1
   DLAMDA( K ) = D( PJ )
   W( K ) = Z( PJ )
   INDXP( K ) = PJ
!
!     Count up the total number of the various types of columns, then
!     form a permutation which positions the four column types into
!     four uniform groups (although one or more of these groups may be
!     empty).
!
   DO 110 J = 1, 4
      CTOT( J ) = 0
110 CONTINUE
   DO 120 J = 1, N
      CT = COLTYP( J )
      CTOT( CT ) = CTOT( CT ) + 1
120 CONTINUE
!
!     PSM(*) = Position in SubMatrix (of types 1 through 4)
!
   PSM( 1 ) = 1
   PSM( 2 ) = 1 + CTOT( 1 )
   PSM( 3 ) = PSM( 2 ) + CTOT( 2 )
   PSM( 4 ) = PSM( 3 ) + CTOT( 3 )
   K = N - CTOT( 4 )
!
!     Fill out the INDXC array so that the permutation which it induces
!     will place all type-1 columns first, all type-2 columns next,
!     then all type-3's, and finally all type-4's.
!
   DO 130 J = 1, N
      JS = INDXP( J )
      CT = COLTYP( JS )
      INDX( PSM( CT ) ) = JS
      INDXC( PSM( CT ) ) = J
      PSM( CT ) = PSM( CT ) + 1
130 CONTINUE
!
!     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
!     and Q2 respectively.  The eigenvalues/vectors which were not
!     deflated go into the first K slots of DLAMDA and Q2 respectively,
!     while those which were deflated go into the last N - K slots.
!
   I = 1
   IQ1 = 1
   IQ2 = 1 + ( CTOT( 1 )+CTOT( 2 ) )*N1
   DO 140 J = 1, CTOT( 1 )
      JS = INDX( I )
      CALL DCOPY( N1, Q( 1, JS ), 1, Q2( IQ1 ), 1 )
      Z( I ) = D( JS )
      I = I + 1
      IQ1 = IQ1 + N1
140 CONTINUE
!
   DO 150 J = 1, CTOT( 2 )
      JS = INDX( I )
      CALL DCOPY( N1, Q( 1, JS ), 1, Q2( IQ1 ), 1 )
      CALL DCOPY( N2, Q( N1+1, JS ), 1, Q2( IQ2 ), 1 )
      Z( I ) = D( JS )
      I = I + 1
      IQ1 = IQ1 + N1
      IQ2 = IQ2 + N2
150 CONTINUE
!
   DO 160 J = 1, CTOT( 3 )
      JS = INDX( I )
      CALL DCOPY( N2, Q( N1+1, JS ), 1, Q2( IQ2 ), 1 )
      Z( I ) = D( JS )
      I = I + 1
      IQ2 = IQ2 + N2
160 CONTINUE
!
   IQ1 = IQ2
   DO 170 J = 1, CTOT( 4 )
      JS = INDX( I )
      CALL DCOPY( N, Q( 1, JS ), 1, Q2( IQ2 ), 1 )
      IQ2 = IQ2 + N
      Z( I ) = D( JS )
      I = I + 1
170 CONTINUE
!
!     The deflated eigenvalues and their corresponding vectors go back
!     into the last N - K slots of D and Q respectively.
!
   IF( K.LT.N ) THEN
      CALL DLACPY( 'A', N, CTOT( 4 ), Q2( IQ1 ), N, &
           Q( 1, K+1 ), LDQ )
      CALL DCOPY( N-K, Z( K+1 ), 1, D( K+1 ), 1 )
   END IF
!
!     Copy CTOT into COLTYP for referencing in DLAED3.
!
   DO 180 J = 1, 4
      COLTYP( J ) = CTOT( J )
180 CONTINUE
!
190 CONTINUE
   RETURN
!
!     End of DLAED2
!
 END SUBROUTINE DLAED2
!
!=
!
!> \brief \b DLAED3 used by sstedc. Finds the roots of the secular equation and updates the eigenvectors. Used when the original matrix is tridiagonal.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED3 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed3.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed3.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed3.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMDA, Q2, INDX,
!                          CTOT, W, S, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDQ, N, N1
!       DOUBLE PRECISION   RHO
!       ..
!       .. Array Arguments ..
!       INTEGER            CTOT( * ), INDX( * )
!       DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),
!      $                   S( * ), W( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED3 finds the roots of the secular equation, as defined by the
!> values in D, W, and RHO, between 1 and K.  It makes the
!> appropriate calls to DLAED4 and then updates the eigenvectors by
!> multiplying the matrix of eigenvectors of the pair of eigensystems
!> being combined by the matrix of eigenvectors of the K-by-K system
!> which is solved here.
!>
!> This code makes very mild assumptions about floating point
!> arithmetic. It will work on machines with a guard digit in
!> add/subtract, or on those binary machines without guard digits
!> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!> It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of terms in the rational function to be solved by
!>          DLAED4.  K >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns in the Q matrix.
!>          N >= K (deflation may result in N>K).
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>          The location of the last eigenvalue in the leading submatrix.
!>          min(1,N) <= N1 <= N/2.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          D(I) contains the updated eigenvalues for
!>          1 <= I <= K.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!>          Initially the first K columns are used as workspace.
!>          On output the columns 1 to K contain
!>          the updated eigenvectors.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>          The value of the parameter in the rank one update equation.
!>          RHO >= 0 required.
!> \endverbatim
!>
!> \param[in,out] DLAMDA
!> \verbatim
!>          DLAMDA is DOUBLE PRECISION array, dimension (K)
!>          The first K elements of this array contain the old roots
!>          of the deflated updating problem.  These are the poles
!>          of the secular equation. May be changed on output by
!>          having lowest order bit set to zero on Cray X-MP, Cray Y-MP,
!>          Cray-2, or Cray C-90, as described above.
!> \endverbatim
!>
!> \param[in] Q2
!> \verbatim
!>          Q2 is DOUBLE PRECISION array, dimension (LDQ2, N)
!>          The first K columns of this matrix contain the non-deflated
!>          eigenvectors for the split problem.
!> \endverbatim
!>
!> \param[in] INDX
!> \verbatim
!>          INDX is INTEGER array, dimension (N)
!>          The permutation used to arrange the columns of the deflated
!>          Q matrix into three groups (see DLAED2).
!>          The rows of the eigenvectors found by DLAED4 must be likewise
!>          permuted before the matrix multiply can take place.
!> \endverbatim
!>
!> \param[in] CTOT
!> \verbatim
!>          CTOT is INTEGER array, dimension (4)
!>          A count of the total number of the various types of columns
!>          in Q, as described in INDX.  The fourth column type is any
!>          column which has been deflated.
!> \endverbatim
!>
!> \param[in,out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (K)
!>          The first K elements of this array contain the components
!>          of the deflation-adjusted updating vector. Destroyed on
!>          output.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (N1 + 1)*K
!>          Will contain the eigenvectors of the repaired matrix which
!>          will be multiplied by the previously accumulated eigenvectors
!>          to update the system.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, an eigenvalue did not converge
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA \n
!>  Modified by Francoise Tisseur, University of Tennessee
!>
!  =====================================================================
 SUBROUTINE DLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMDA, Q2, INDX,&
      CTOT, W, S, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
   INTEGER            INFO, K, LDQ, N, N1
   DOUBLE PRECISION   RHO
!     ..
!     .. Array Arguments ..
   INTEGER            CTOT( * ), INDX( * )
   DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),&
        S( * ), W( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ONE, ZERO
   PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
   INTEGER            I, II, IQ2, J, N12, N2, N23
   DOUBLE PRECISION   TEMP
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMC3, DNRM2
!      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DGEMM, DLACPY, DLAED4, DLASET, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   INFO = 0
!
   IF( K.LT.0 ) THEN
      INFO = -1
   ELSE IF( N.LT.K ) THEN
      INFO = -2
   ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
      INFO = -6
   END IF
   IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'DLAED3', -INFO )
      RETURN
   END IF
!
!     Quick return if possible
!
   IF( K.EQ.0 ) RETURN
!
!     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DLAMDA(I) if it is 1; this makes the subsequent
!     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DLAMDA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DLAMDA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
!     this code.
!
   DO 10 I = 1, K
      DLAMDA( I ) = DLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
10 CONTINUE
!
   DO 20 J = 1, K
      CALL DLAED4( K, J, DLAMDA, W, Q( 1, J ), RHO, D( J ), INFO )
!
!        If the zero finder fails, the computation is terminated.
!
      IF( INFO.NE.0 ) GO TO 120
20 CONTINUE
!
   IF( K.EQ.1 ) GO TO 110
   IF( K.EQ.2 ) THEN
      DO 30 J = 1, K
         W( 1 ) = Q( 1, J )
         W( 2 ) = Q( 2, J )
         II = INDX( 1 )
         Q( 1, J ) = W( II )
         II = INDX( 2 )
         Q( 2, J ) = W( II )
30       CONTINUE
      GO TO 110
   END IF
!
!     Compute updated W.
!
   CALL DCOPY( K, W, 1, S, 1 )
!
!     Initialize W(I) = Q(I,I)
!
   CALL DCOPY( K, Q, LDQ+1, W, 1 )
   DO 60 J = 1, K
      DO 40 I = 1, J - 1
         W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
40    CONTINUE
      DO 50 I = J + 1, K
         W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
50    CONTINUE
60 CONTINUE
   DO 70 I = 1, K
      W( I ) = SIGN( SQRT( -W( I ) ), S( I ) )
70 CONTINUE
!
!     Compute eigenvectors of the modified rank-1 modification.
!
   DO 100 J = 1, K
      DO 80 I = 1, K
         S( I ) = W( I ) / Q( I, J )
80    CONTINUE
      TEMP = DNRM2( K, S, 1 )
      DO 90 I = 1, K
         II = INDX( I )
         Q( I, J ) = S( II ) / TEMP
90    CONTINUE
100 CONTINUE
!
!     Compute the updated eigenvectors.
!
110 CONTINUE
!
   N2 = N - N1
   N12 = CTOT( 1 ) + CTOT( 2 )
   N23 = CTOT( 2 ) + CTOT( 3 )
!
   CALL DLACPY( 'A', N23, K, Q( CTOT( 1 )+1, 1 ), LDQ, S, N23 )
   IQ2 = N1*N12 + 1
   IF( N23.NE.0 ) THEN
      CALL DGEMM( 'N', 'N', N2, K, N23, ONE, Q2( IQ2 ), N2, S, N23,&
           ZERO, Q( N1+1, 1 ), LDQ )
   ELSE
      CALL DLASET( 'A', N2, K, ZERO, ZERO, Q( N1+1, 1 ), LDQ )
   END IF
!
   CALL DLACPY( 'A', N12, K, Q, LDQ, S, N12 )
   IF( N12.NE.0 ) THEN
      CALL DGEMM( 'N', 'N', N1, K, N12, ONE, Q2, N1, S, N12, ZERO, Q,&
           LDQ )
   ELSE
      CALL DLASET( 'A', N1, K, ZERO, ZERO, Q( 1, 1 ), LDQ )
   END IF
!
!
120 CONTINUE
   RETURN
!
!     End of DLAED3
!
 END SUBROUTINE DLAED3
!
!=
!
!> \brief \b DLAED4 used by sstedc. Finds a single root of the secular equation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED4 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed4.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed4.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed4.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            I, INFO, N
!       DOUBLE PRECISION   DLAM, RHO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), DELTA( * ), Z( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This subroutine computes the I-th updated eigenvalue of a symmetric
!> rank-one modification to a diagonal matrix whose elements are
!> given in the array d, and that
!>
!>            D(i) < D(j)  for  i < j
!>
!> and that RHO > 0.  This is arranged by the calling routine, and is
!> no loss in generality.  The rank-one modified system is thus
!>
!>            diag( D )  +  RHO * Z * Z_transpose.
!>
!> where we assume the Euclidean norm of Z is 1.
!>
!> The method consists of approximating the rational functions in the
!> secular equation by simpler interpolating rational functions.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The length of all arrays.
!> \endverbatim
!>
!> \param[in] I
!> \verbatim
!>          I is INTEGER
!>         The index of the eigenvalue to be computed.  1 <= I <= N.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>         The original eigenvalues.  It is assumed that they are in
!>         order, D(I) < D(J)  for I < J.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (N)
!>         The components of the updating vector.
!> \endverbatim
!>
!> \param[out] DELTA
!> \verbatim
!>          DELTA is DOUBLE PRECISION array, dimension (N)
!>         If N .GT. 2, DELTA contains (D(j) - lambda_I) in its  j-th
!>         component.  If N = 1, then DELTA(1) = 1. If N = 2, see DLAED5
!>         for detail. The vector DELTA contains the information necessary
!>         to construct the eigenvectors by DLAED3 and DLAED9.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>         The scalar in the symmetric updating formula.
!> \endverbatim
!>
!> \param[out] DLAM
!> \verbatim
!>          DLAM is DOUBLE PRECISION
!>         The computed lambda_I, the I-th updated eigenvalue.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>         = 0:  successful exit
!>         > 0:  if INFO = 1, the updating process failed.
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  Logical variable ORGATI (origin-at-i?) is used for distinguishing
!>  whether D(i) or D(i+1) is treated as the origin.
!>
!>            ORGATI = .true.    origin at i
!>            ORGATI = .false.   origin at i+1
!>
!>   Logical variable SWTCH3 (switch-for-3-poles?) is for noting
!>   if we are working with THREE poles!
!>
!>   MAXIT is the maximum number of iterations allowed for each
!>   eigenvalue.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>     Ren-Cang Li, Computer Science Division, University of California
!>     at Berkeley, USA
!>
!  =====================================================================
 SUBROUTINE DLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
   INTEGER            I, INFO, N
   DOUBLE PRECISION   DLAM, RHO
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   D( * ), DELTA( * ), Z( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXIT
   PARAMETER          ( MAXIT = 30 )
   DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR, EIGHT, TEN
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,&
        THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0,&
        TEN = 10.0D0 )
!     ..
!     .. Local Scalars ..
   LOGICAL            ORGATI, SWTCH, SWTCH3
   INTEGER            II, IIM1, IIP1, IP1, ITER, J, NITER
   DOUBLE PRECISION   A, B, C, DEL, DLTLB, DLTUB, DPHI, DPSI, DW,&
        EPS, ERRETM, ETA, MIDPT, PHI, PREW, PSI,&
        RHOINV, TAU, TEMP, TEMP1, W
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   ZZ( 3 )
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLAED5, DLAED6
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Since this routine is called in an inner loop, we do no argument
!     checking.
!
!     Quick return for N=1 and 2.
!
   INFO = 0
   IF( N.EQ.1 ) THEN
!
!         Presumably, I=1 upon entry
!
      DLAM = D( 1 ) + RHO*Z( 1 )*Z( 1 )
      DELTA( 1 ) = ONE
      RETURN
   END IF
   IF( N.EQ.2 ) THEN
      CALL DLAED5( I, D, Z, DELTA, RHO, DLAM )
      RETURN
   END IF
!
!     Compute machine epsilon
!
   EPS = DLAMCH( 'Epsilon' )
   RHOINV = ONE / RHO
!
!     The case I = N
!
   IF( I.EQ.N ) THEN
!
!        Initialize some basic variables
!
      II = N - 1
      NITER = 1
!
!        Calculate initial guess
!
      MIDPT = RHO / TWO
!
!        If ||Z||_2 is not one, then TEMP should be set to
!        RHO * ||Z||_2^2 / TWO
!
      DO 10 J = 1, N
         DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
10    CONTINUE
!
      PSI = ZERO
      DO 20 J = 1, N - 2
         PSI = PSI + Z( J )*Z( J ) / DELTA( J )
20    CONTINUE
!
      C = RHOINV + PSI
      W = C + Z( II )*Z( II ) / DELTA( II ) +&
           Z( N )*Z( N ) / DELTA( N )
!
      IF( W.LE.ZERO ) THEN
         TEMP = Z( N-1 )*Z( N-1 ) / ( D( N )-D( N-1 )+RHO ) +&
              Z( N )*Z( N ) / RHO
         IF( C.LE.TEMP ) THEN
            TAU = RHO
         ELSE
            DEL = D( N ) - D( N-1 )
            A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
            B = Z( N )*Z( N )*DEL
            IF( A.LT.ZERO ) THEN
               TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
            ELSE
               TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
            END IF
         END IF
!
!           It can be proved that
!               D(N)+RHO/2 <= LAMBDA(N) < D(N)+TAU <= D(N)+RHO
!
         DLTLB = MIDPT
         DLTUB = RHO
      ELSE
         DEL = D( N ) - D( N-1 )
         A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
         B = Z( N )*Z( N )*DEL
         IF( A.LT.ZERO ) THEN
            TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
         ELSE
            TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
         END IF
!
!           It can be proved that
!               D(N) < D(N)+TAU < LAMBDA(N) < D(N)+RHO/2
!
         DLTLB = ZERO
         DLTUB = MIDPT
      END IF
!
      DO 30 J = 1, N
         DELTA( J ) = ( D( J )-D( I ) ) - TAU
30    CONTINUE
!
!        Evaluate PSI and the derivative DPSI
!
      DPSI = ZERO
      PSI = ZERO
      ERRETM = ZERO
      DO 40 J = 1, II
         TEMP = Z( J ) / DELTA( J )
         PSI = PSI + Z( J )*TEMP
         DPSI = DPSI + TEMP*TEMP
         ERRETM = ERRETM + PSI
40    CONTINUE
      ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
      TEMP = Z( N ) / DELTA( N )
      PHI = Z( N )*TEMP
      DPHI = TEMP*TEMP
      ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +&
           ABS( TAU )*( DPSI+DPHI )
!
      W = RHOINV + PHI + PSI
!
!        Test for convergence
!
      IF( ABS( W ).LE.EPS*ERRETM ) THEN
         DLAM = D( I ) + TAU
         GO TO 250
      END IF
!
      IF( W.LE.ZERO ) THEN
         DLTLB = MAX( DLTLB, TAU )
      ELSE
         DLTUB = MIN( DLTUB, TAU )
      END IF
!
!        Calculate the new step
!
      NITER = NITER + 1
      C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
      A = ( DELTA( N-1 )+DELTA( N ) )*W -&
           DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
      B = DELTA( N-1 )*DELTA( N )*W
      IF( C.LT.ZERO ) C = ABS( C )
      IF( C.EQ.ZERO ) THEN
!          ETA = B/A
!           ETA = RHO - TAU
         ETA = DLTUB - TAU
      ELSE IF( A.GE.ZERO ) THEN
         ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
      ELSE
         ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
      END IF
!
!        Note, eta should be positive if w is negative, and
!        eta should be negative otherwise. However,
!        if for some reason caused by roundoff, eta*w > 0,
!        we simply use one Newton step instead. This way
!        will guarantee eta*w < 0.
!
      IF( W*ETA.GT.ZERO ) ETA = -W / ( DPSI+DPHI )
      TEMP = TAU + ETA
      IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
         IF( W.LT.ZERO ) THEN
            ETA = ( DLTUB-TAU ) / TWO
         ELSE
            ETA = ( DLTLB-TAU ) / TWO
         END IF
      END IF
      DO 50 J = 1, N
         DELTA( J ) = DELTA( J ) - ETA
50    CONTINUE
!
      TAU = TAU + ETA
!
!        Evaluate PSI and the derivative DPSI
!
      DPSI = ZERO
      PSI = ZERO
      ERRETM = ZERO
      DO 60 J = 1, II
         TEMP = Z( J ) / DELTA( J )
         PSI = PSI + Z( J )*TEMP
         DPSI = DPSI + TEMP*TEMP
         ERRETM = ERRETM + PSI
60    CONTINUE
      ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
      TEMP = Z( N ) / DELTA( N )
      PHI = Z( N )*TEMP
      DPHI = TEMP*TEMP
      ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +&
           ABS( TAU )*( DPSI+DPHI )
!
      W = RHOINV + PHI + PSI
!
!        Main loop to update the values of the array   DELTA
!
      ITER = NITER + 1
!
      DO 90 NITER = ITER, MAXIT
!
!           Test for convergence
!
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            DLAM = D( I ) + TAU
            GO TO 250
         END IF
!
         IF( W.LE.ZERO ) THEN
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         END IF
!
!           Calculate the new step
!
         C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
         A = ( DELTA( N-1 )+DELTA( N ) )*W -&
              DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
         B = DELTA( N-1 )*DELTA( N )*W
         IF( A.GE.ZERO ) THEN
            ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
!
!           Note, eta should be positive if w is negative, and
!           eta should be negative otherwise. However,
!           if for some reason caused by roundoff, eta*w > 0,
!           we simply use one Newton step instead. This way
!           will guarantee eta*w < 0.
!
         IF( W*ETA.GT.ZERO )&
              ETA = -W / ( DPSI+DPHI )
         TEMP = TAU + ETA
         IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            END IF
         END IF
         DO 70 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
70       CONTINUE
!
         TAU = TAU + ETA
!
!           Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 80 J = 1, II
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
80       CONTINUE
         ERRETM = ABS( ERRETM )
!
!           Evaluate PHI and the derivative DPHI
!
         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +&
              ABS( TAU )*( DPSI+DPHI )
!
         W = RHOINV + PHI + PSI
90    CONTINUE
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
      INFO = 1
      DLAM = D( I ) + TAU
      GO TO 250
!
!        End for the case I = N
!
   ELSE
!
!        The case for I < N
!
      NITER = 1
      IP1 = I + 1
!
!        Calculate initial guess
!
      DEL = D( IP1 ) - D( I )
      MIDPT = DEL / TWO
      DO 100 J = 1, N
         DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
100   CONTINUE
!
      PSI = ZERO
      DO 110 J = 1, I - 1
         PSI = PSI + Z( J )*Z( J ) / DELTA( J )
110   CONTINUE
!
      PHI = ZERO
      DO 120 J = N, I + 2, -1
         PHI = PHI + Z( J )*Z( J ) / DELTA( J )
120   CONTINUE
      C = RHOINV + PSI + PHI
      W = C + Z( I )*Z( I ) / DELTA( I ) +&
           Z( IP1 )*Z( IP1 ) / DELTA( IP1 )
!
      IF( W.GT.ZERO ) THEN
!
!           d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
!
!           We choose d(i) as origin.
!
         ORGATI = .TRUE.
         A = C*DEL + Z( I )*Z( I ) + Z( IP1 )*Z( IP1 )
         B = Z( I )*Z( I )*DEL
         IF( A.GT.ZERO ) THEN
            TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
         ELSE
            TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         END IF
         DLTLB = ZERO
         DLTUB = MIDPT
      ELSE
!
!           (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
!
!           We choose d(i+1) as origin.
!
         ORGATI = .FALSE.
         A = C*DEL - Z( I )*Z( I ) - Z( IP1 )*Z( IP1 )
         B = Z( IP1 )*Z( IP1 )*DEL
         IF( A.LT.ZERO ) THEN
            TAU = TWO*B / ( A-SQRT( ABS( A*A+FOUR*B*C ) ) )
         ELSE
            TAU = -( A+SQRT( ABS( A*A+FOUR*B*C ) ) ) / ( TWO*C )
         END IF
         DLTLB = -MIDPT
         DLTUB = ZERO
      END IF
!
      IF( ORGATI ) THEN
         DO 130 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - TAU
130      CONTINUE
      ELSE
         DO 140 J = 1, N
            DELTA( J ) = ( D( J )-D( IP1 ) ) - TAU
140      CONTINUE
      END IF
      IF( ORGATI ) THEN
         II = I
      ELSE
         II = I + 1
      END IF
      IIM1 = II - 1
      IIP1 = II + 1
!
!        Evaluate PSI and the derivative DPSI
!
      DPSI = ZERO
      PSI = ZERO
      ERRETM = ZERO
      DO 150 J = 1, IIM1
         TEMP = Z( J ) / DELTA( J )
         PSI = PSI + Z( J )*TEMP
         DPSI = DPSI + TEMP*TEMP
         ERRETM = ERRETM + PSI
150   CONTINUE
      ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
      DPHI = ZERO
      PHI = ZERO
      DO 160 J = N, IIP1, -1
         TEMP = Z( J ) / DELTA( J )
         PHI = PHI + Z( J )*TEMP
         DPHI = DPHI + TEMP*TEMP
         ERRETM = ERRETM + PHI
160   CONTINUE
!
      W = RHOINV + PHI + PSI
!
!        W is the value of the secular function with
!        its ii-th element removed.
!
      SWTCH3 = .FALSE.
      IF( ORGATI ) THEN
         IF( W.LT.ZERO ) SWTCH3 = .TRUE.
      ELSE
         IF( W.GT.ZERO ) SWTCH3 = .TRUE.
      END IF
      IF( II.EQ.1 .OR. II.EQ.N ) SWTCH3 = .FALSE.
!
      TEMP = Z( II ) / DELTA( II )
      DW = DPSI + DPHI + TEMP*TEMP
      TEMP = Z( II )*TEMP
      W = W + TEMP
      ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +&
           THREE*ABS( TEMP ) + ABS( TAU )*DW
!
!        Test for convergence
!
      IF( ABS( W ).LE.EPS*ERRETM ) THEN
         IF( ORGATI ) THEN
            DLAM = D( I ) + TAU
         ELSE
            DLAM = D( IP1 ) + TAU
         END IF
         GO TO 250
      END IF
!
      IF( W.LE.ZERO ) THEN
         DLTLB = MAX( DLTLB, TAU )
      ELSE
         DLTUB = MIN( DLTUB, TAU )
      END IF
!
!        Calculate the new step
!
      NITER = NITER + 1
      IF( .NOT.SWTCH3 ) THEN
         IF( ORGATI ) THEN
            C = W - DELTA( IP1 )*DW - ( D( I )-D( IP1 ) )*&
                 ( Z( I ) / DELTA( I ) )**2
         ELSE
            C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )*&
                 ( Z( IP1 ) / DELTA( IP1 ) )**2
         END IF
         A = ( DELTA( I )+DELTA( IP1 ) )*W -&
              DELTA( I )*DELTA( IP1 )*DW
         B = DELTA( I )*DELTA( IP1 )*W
         IF( C.EQ.ZERO ) THEN
            IF( A.EQ.ZERO ) THEN
               IF( ORGATI ) THEN
                  A = Z( I )*Z( I ) + DELTA( IP1 )*DELTA( IP1 )*&
                       ( DPSI+DPHI )
               ELSE
                  A = Z( IP1 )*Z( IP1 ) + DELTA( I )*DELTA( I )*&
                       ( DPSI+DPHI )
               END IF
            END IF
            ETA = B / A
         ELSE IF( A.LE.ZERO ) THEN
            ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
      ELSE
!
!           Interpolation using THREE most relevant poles
!
         TEMP = RHOINV + PSI + PHI
         IF( ORGATI ) THEN
            TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
            TEMP1 = TEMP1*TEMP1
            C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) -&
                 ( D( IIM1 )-D( IIP1 ) )*TEMP1
            ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
            ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*&
                 ( ( DPSI-TEMP1 )+DPHI )
         ELSE
            TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
            TEMP1 = TEMP1*TEMP1
            C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) -&
                 ( D( IIP1 )-D( IIM1 ) )*TEMP1
            ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*&
                 ( DPSI+( DPHI-TEMP1 ) )
            ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
         END IF
         ZZ( 2 ) = Z( II )*Z( II )
         CALL DLAED6( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA,&
              INFO )
         IF( INFO.NE.0 ) GO TO 250
      END IF
!
!        Note, eta should be positive if w is negative, and
!        eta should be negative otherwise. However,
!        if for some reason caused by roundoff, eta*w > 0,
!        we simply use one Newton step instead. This way
!        will guarantee eta*w < 0.
!
      IF( W*ETA.GE.ZERO ) ETA = -W / DW
      TEMP = TAU + ETA
      IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
         IF( W.LT.ZERO ) THEN
            ETA = ( DLTUB-TAU ) / TWO
         ELSE
            ETA = ( DLTLB-TAU ) / TWO
         END IF
      END IF
!
      PREW = W
!
      DO 180 J = 1, N
         DELTA( J ) = DELTA( J ) - ETA
180   CONTINUE
!
!        Evaluate PSI and the derivative DPSI
!
      DPSI = ZERO
      PSI = ZERO
      ERRETM = ZERO
      DO 190 J = 1, IIM1
         TEMP = Z( J ) / DELTA( J )
         PSI = PSI + Z( J )*TEMP
         DPSI = DPSI + TEMP*TEMP
         ERRETM = ERRETM + PSI
190   CONTINUE
      ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
      DPHI = ZERO
      PHI = ZERO
      DO 200 J = N, IIP1, -1
         TEMP = Z( J ) / DELTA( J )
         PHI = PHI + Z( J )*TEMP
         DPHI = DPHI + TEMP*TEMP
         ERRETM = ERRETM + PHI
200   CONTINUE
!
      TEMP = Z( II ) / DELTA( II )
      DW = DPSI + DPHI + TEMP*TEMP
      TEMP = Z( II )*TEMP
      W = RHOINV + PHI + PSI + TEMP
      ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +&
           THREE*ABS( TEMP ) + ABS( TAU+ETA )*DW
!
      SWTCH = .FALSE.
      IF( ORGATI ) THEN
         IF( -W.GT.ABS( PREW ) / TEN ) SWTCH = .TRUE.
      ELSE
         IF( W.GT.ABS( PREW ) / TEN )  SWTCH = .TRUE.
      END IF
!
      TAU = TAU + ETA
!
!        Main loop to update the values of the array   DELTA
!
      ITER = NITER + 1
!
      DO 240 NITER = ITER, MAXIT
!
!           Test for convergence
!
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            IF( ORGATI ) THEN
               DLAM = D( I ) + TAU
            ELSE
               DLAM = D( IP1 ) + TAU
            END IF
            GO TO 250
         END IF
!
         IF( W.LE.ZERO ) THEN
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         END IF
!
!           Calculate the new step
!
         IF( .NOT.SWTCH3 ) THEN
            IF( .NOT.SWTCH ) THEN
               IF( ORGATI ) THEN
                  C = W - DELTA( IP1 )*DW -&
                       ( D( I )-D( IP1 ) )*( Z( I ) / DELTA( I ) )**2
               ELSE
                  C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )*&
                       ( Z( IP1 ) / DELTA( IP1 ) )**2
               END IF
            ELSE
               TEMP = Z( II ) / DELTA( II )
               IF( ORGATI ) THEN
                  DPSI = DPSI + TEMP*TEMP
               ELSE
                  DPHI = DPHI + TEMP*TEMP
               END IF
               C = W - DELTA( I )*DPSI - DELTA( IP1 )*DPHI
            END IF
            A = ( DELTA( I )+DELTA( IP1 ) )*W -&
                 DELTA( I )*DELTA( IP1 )*DW
            B = DELTA( I )*DELTA( IP1 )*W
            IF( C.EQ.ZERO ) THEN
               IF( A.EQ.ZERO ) THEN
                  IF( .NOT.SWTCH ) THEN
                     IF( ORGATI ) THEN
                        A = Z( I )*Z( I ) + DELTA( IP1 )*&
                             DELTA( IP1 )*( DPSI+DPHI )
                     ELSE
                        A = Z( IP1 )*Z( IP1 ) +&
                             DELTA( I )*DELTA( I )*( DPSI+DPHI )
                     END IF
                  ELSE
                     A = DELTA( I )*DELTA( I )*DPSI +&
                          DELTA( IP1 )*DELTA( IP1 )*DPHI
                  END IF
               END IF
               ETA = B / A
            ELSE IF( A.LE.ZERO ) THEN
               ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
         ELSE
!
!              Interpolation using THREE most relevant poles
!
            TEMP = RHOINV + PSI + PHI
            IF( SWTCH ) THEN
               C = TEMP - DELTA( IIM1 )*DPSI - DELTA( IIP1 )*DPHI
               ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*DPSI
               ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*DPHI
            ELSE
               IF( ORGATI ) THEN
                  TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
                  TEMP1 = TEMP1*TEMP1
                  C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) -&
                       ( D( IIM1 )-D( IIP1 ) )*TEMP1
                  ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
                  ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*&
                       ( ( DPSI-TEMP1 )+DPHI )
               ELSE
                  TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
                  TEMP1 = TEMP1*TEMP1
                  C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) -&
                       ( D( IIP1 )-D( IIM1 ) )*TEMP1
                  ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*&
                       ( DPSI+( DPHI-TEMP1 ) )
                  ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
               END IF
            END IF
            CALL DLAED6( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA,&
                 INFO )
            IF( INFO.NE.0 ) GO TO 250
         END IF
!
!           Note, eta should be positive if w is negative, and
!           eta should be negative otherwise. However,
!           if for some reason caused by roundoff, eta*w > 0,
!           we simply use one Newton step instead. This way
!           will guarantee eta*w < 0.
!
         IF( W*ETA.GE.ZERO ) ETA = -W / DW
         TEMP = TAU + ETA
         IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            END IF
         END IF
!
         DO 210 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
210      CONTINUE
!
         TAU = TAU + ETA
         PREW = W
!
!           Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 220 J = 1, IIM1
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
220      CONTINUE
         ERRETM = ABS( ERRETM )
!
!           Evaluate PHI and the derivative DPHI
!
         DPHI = ZERO
         PHI = ZERO
         DO 230 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
230      CONTINUE
!
         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = RHOINV + PHI + PSI + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +&
              THREE*ABS( TEMP ) + ABS( TAU )*DW
         IF( W*PREW.GT.ZERO .AND. ABS( W ).GT.ABS( PREW ) / TEN ) &
              SWTCH = .NOT.SWTCH
!
240   CONTINUE
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
      INFO = 1
      IF( ORGATI ) THEN
         DLAM = D( I ) + TAU
      ELSE
         DLAM = D( IP1 ) + TAU
      END IF
!
   END IF
!
250 CONTINUE
!
   RETURN
!
!     End of DLAED4
!
 END SUBROUTINE DLAED4
!
!=
!
!> \brief \b DLAED5 used by sstedc. Solves the 2-by-2 secular equation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED5 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed5.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed5.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed5.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED5( I, D, Z, DELTA, RHO, DLAM )
! 
!       .. Scalar Arguments ..
!       INTEGER            I
!       DOUBLE PRECISION   DLAM, RHO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( 2 ), DELTA( 2 ), Z( 2 )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This subroutine computes the I-th eigenvalue of a symmetric rank-one
!> modification of a 2-by-2 diagonal matrix
!>
!>            diag( D )  +  RHO * Z * transpose(Z) .
!>
!> The diagonal elements in the array D are assumed to satisfy
!>
!>            D(i) < D(j)  for  i < j .
!>
!> We also assume RHO > 0 and that the Euclidean norm of the vector
!> Z is one.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] I
!> \verbatim
!>          I is INTEGER
!>         The index of the eigenvalue to be computed.  I = 1 or I = 2.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (2)
!>         The original eigenvalues.  We assume D(1) < D(2).
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (2)
!>         The components of the updating vector.
!> \endverbatim
!>
!> \param[out] DELTA
!> \verbatim
!>          DELTA is DOUBLE PRECISION array, dimension (2)
!>         The vector DELTA contains the information necessary
!>         to construct the eigenvectors.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>         The scalar in the symmetric updating formula.
!> \endverbatim
!>
!> \param[out] DLAM
!> \verbatim
!>          DLAM is DOUBLE PRECISION
!>         The computed lambda_I, the I-th updated eigenvalue.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>     Ren-Cang Li, Computer Science Division, University of California
!>     at Berkeley, USA
!>
!  =====================================================================
 SUBROUTINE DLAED5( I, D, Z, DELTA, RHO, DLAM )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
   INTEGER            I
   DOUBLE PRECISION   DLAM, RHO
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   D( 2 ), DELTA( 2 ), Z( 2 )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION   ZERO, ONE, TWO, FOUR
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,&
        FOUR = 4.0D0 )
!     ..
!     .. Local Scalars ..
   DOUBLE PRECISION   B, C, DEL, TAU, TEMP, W
!     ..
!     .. Intrinsic Functions ..
!   INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
   DEL = D( 2 ) - D( 1 )
   IF( I.EQ.1 ) THEN
      W = ONE + TWO*RHO*( Z( 2 )*Z( 2 )-Z( 1 )*Z( 1 ) ) / DEL
      IF( W.GT.ZERO ) THEN
         B = DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
         C = RHO*Z( 1 )*Z( 1 )*DEL
!
!           B > ZERO, always
!
         TAU = TWO*C / ( B+SQRT( ABS( B*B-FOUR*C ) ) )
         DLAM = D( 1 ) + TAU
         DELTA( 1 ) = -Z( 1 ) / TAU
         DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
      ELSE
         B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
         C = RHO*Z( 2 )*Z( 2 )*DEL
         IF( B.GT.ZERO ) THEN
            TAU = -TWO*C / ( B+SQRT( B*B+FOUR*C ) )
         ELSE
            TAU = ( B-SQRT( B*B+FOUR*C ) ) / TWO
         END IF
         DLAM = D( 2 ) + TAU
         DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
         DELTA( 2 ) = -Z( 2 ) / TAU
      END IF
      TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
      DELTA( 1 ) = DELTA( 1 ) / TEMP
      DELTA( 2 ) = DELTA( 2 ) / TEMP
   ELSE
!
!     Now I=2
!
      B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
      C = RHO*Z( 2 )*Z( 2 )*DEL
      IF( B.GT.ZERO ) THEN
         TAU = ( B+SQRT( B*B+FOUR*C ) ) / TWO
      ELSE
         TAU = TWO*C / ( -B+SQRT( B*B+FOUR*C ) )
      END IF
      DLAM = D( 2 ) + TAU
      DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
      DELTA( 2 ) = -Z( 2 ) / TAU
      TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
      DELTA( 1 ) = DELTA( 1 ) / TEMP
      DELTA( 2 ) = DELTA( 2 ) / TEMP
   END IF
   RETURN
!
!     End OF DLAED5
!
 END SUBROUTINE DLAED5
!=
!> \brief \b DLAED6 used by sstedc. Computes one Newton step in solution of the secular equation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED6 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed6.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed6.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed6.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO )
! 
!       .. Scalar Arguments ..
!       LOGICAL            ORGATI
!       INTEGER            INFO, KNITER
!       DOUBLE PRECISION   FINIT, RHO, TAU
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( 3 ), Z( 3 )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED6 computes the positive or negative root (closest to the origin)
!> of
!>                  z(1)        z(2)        z(3)
!> f(x) =   rho + --------- + ---------- + ---------
!>                 d(1)-x      d(2)-x      d(3)-x
!>
!> It is assumed that
!>
!>       if ORGATI = .true. the root is between d(2) and d(3);
!>       otherwise it is between d(1) and d(2)
!>
!> This routine will be called by DLAED4 when necessary. In most cases,
!> the root sought is the smallest in magnitude, though it might not be
!> in some extremely rare situations.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] KNITER
!> \verbatim
!>          KNITER is INTEGER
!>               Refer to DLAED4 for its significance.
!> \endverbatim
!>
!> \param[in] ORGATI
!> \verbatim
!>          ORGATI is LOGICAL
!>               If ORGATI is true, the needed root is between d(2) and
!>               d(3); otherwise it is between d(1) and d(2).  See
!>               DLAED4 for further details.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>               Refer to the equation f(x) above.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (3)
!>               D satisfies d(1) < d(2) < d(3).
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (3)
!>               Each of the elements in z must be positive.
!> \endverbatim
!>
!> \param[in] FINIT
!> \verbatim
!>          FINIT is DOUBLE PRECISION
!>               The value of f at 0. It is more accurate than the one
!>               evaluated inside this routine (if someone wants to do
!>               so).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION
!>               The root of the equation f(x).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>               = 0: successful exit
!>               > 0: if INFO = 1, failure to converge
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2015
!
!> \ingroup auxOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  10/02/03: This version has a few statements commented out for thread
!>  safety (machine parameters are computed on each entry). SJH.
!>
!>  05/10/06: Modified from a new version of Ren-Cang Li, use
!>     Gragg-Thornton-Warner cubic convergent scheme for better stability.
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>     Ren-Cang Li, Computer Science Division, University of California
!>     at Berkeley, USA
!>
!  =====================================================================
 SUBROUTINE DLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO )
!
!  -- LAPACK computational routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
   LOGICAL            ORGATI
   INTEGER            INFO, KNITER
   DOUBLE PRECISION   FINIT, RHO, TAU
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION   D( 3 ), Z( 3 )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   INTEGER            MAXIT
   PARAMETER          ( MAXIT = 40 )
   DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR, EIGHT
   PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,&
        THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0 )
!     ..
!     .. External Functions ..
!   DOUBLE PRECISION   DLAMCH
!   EXTERNAL           DLAMCH
!     ..
!     .. Local Arrays ..
   DOUBLE PRECISION   DSCALE( 3 ), ZSCALE( 3 )
!     ..
!     .. Local Scalars ..
   LOGICAL            SCALE
   INTEGER            I, ITER, NITER
   DOUBLE PRECISION   A, B, BASE, C, DDF, DF, EPS, ERRETM, ETA, F,&
        FC, SCLFAC, SCLINV, SMALL1, SMALL2, SMINV1,&
        SMINV2, TEMP, TEMP1, TEMP2, TEMP3, TEMP4, &
        LBD, UBD
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
   INFO = 0
!
   IF( ORGATI ) THEN
      LBD = D(2)
      UBD = D(3)
   ELSE
      LBD = D(1)
      UBD = D(2)
   END IF
   IF( FINIT .LT. ZERO )THEN
      LBD = ZERO
   ELSE
      UBD = ZERO 
   END IF
!
   NITER = 1
   TAU = ZERO
   IF( KNITER.EQ.2 ) THEN
      IF( ORGATI ) THEN
         TEMP = ( D( 3 )-D( 2 ) ) / TWO
         C = RHO + Z( 1 ) / ( ( D( 1 )-D( 2 ) )-TEMP )
         A = C*( D( 2 )+D( 3 ) ) + Z( 2 ) + Z( 3 )
         B = C*D( 2 )*D( 3 ) + Z( 2 )*D( 3 ) + Z( 3 )*D( 2 )
      ELSE
         TEMP = ( D( 1 )-D( 2 ) ) / TWO
         C = RHO + Z( 3 ) / ( ( D( 3 )-D( 2 ) )-TEMP )
         A = C*( D( 1 )+D( 2 ) ) + Z( 1 ) + Z( 2 )
         B = C*D( 1 )*D( 2 ) + Z( 1 )*D( 2 ) + Z( 2 )*D( 1 )
      END IF
      TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
      A = A / TEMP
      B = B / TEMP
      C = C / TEMP
      IF( C.EQ.ZERO ) THEN
         TAU = B / A
      ELSE IF( A.LE.ZERO ) THEN
         TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
      ELSE
         TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
      END IF
      IF( TAU .LT. LBD .OR. TAU .GT. UBD ) TAU = ( LBD+UBD )/TWO
      IF( D(1).EQ.TAU .OR. D(2).EQ.TAU .OR. D(3).EQ.TAU ) THEN
         TAU = ZERO
      ELSE
         TEMP = FINIT + TAU*Z(1)/( D(1)*( D( 1 )-TAU ) ) +&
              TAU*Z(2)/( D(2)*( D( 2 )-TAU ) ) +&
              TAU*Z(3)/( D(3)*( D( 3 )-TAU ) )
         IF( TEMP .LE. ZERO )THEN
            LBD = TAU
         ELSE
            UBD = TAU
         END IF
         IF( ABS( FINIT ).LE.ABS( TEMP ) ) TAU = ZERO
      END IF
   END IF
!
!     get machine parameters for possible scaling to avoid overflow
!
!     modified by Sven: parameters SMALL1, SMINV1, SMALL2,
!     SMINV2, EPS are not SAVEd anymore between one call to the
!     others but recomputed at each call
!
   EPS = DLAMCH( 'Epsilon' )
   BASE = DLAMCH( 'Base' )
   SMALL1 = BASE**( INT( LOG( DLAMCH( 'SafMin' ) ) / LOG( BASE ) /&
        THREE ) )
   SMINV1 = ONE / SMALL1
   SMALL2 = SMALL1*SMALL1
   SMINV2 = SMINV1*SMINV1
!
!     Determine if scaling of inputs necessary to avoid overflow
!     when computing 1/TEMP**3
!
   IF( ORGATI ) THEN
      TEMP = MIN( ABS( D( 2 )-TAU ), ABS( D( 3 )-TAU ) )
   ELSE
      TEMP = MIN( ABS( D( 1 )-TAU ), ABS( D( 2 )-TAU ) )
   END IF
   SCALE = .FALSE.
   IF( TEMP.LE.SMALL1 ) THEN
      SCALE = .TRUE.
      IF( TEMP.LE.SMALL2 ) THEN
!
!        Scale up by power of radix nearest 1/SAFMIN**(2/3)
!
         SCLFAC = SMINV2
         SCLINV = SMALL2
      ELSE
!
!        Scale up by power of radix nearest 1/SAFMIN**(1/3)
!
         SCLFAC = SMINV1
         SCLINV = SMALL1
      END IF
!
!        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1)
!
      DO 10 I = 1, 3
         DSCALE( I ) = D( I )*SCLFAC
         ZSCALE( I ) = Z( I )*SCLFAC
10    CONTINUE
      TAU = TAU*SCLFAC
      LBD = LBD*SCLFAC
      UBD = UBD*SCLFAC
   ELSE
!
!        Copy D and Z to DSCALE and ZSCALE
!
      DO 20 I = 1, 3
         DSCALE( I ) = D( I )
         ZSCALE( I ) = Z( I )
20    CONTINUE
      END IF
!
      FC = ZERO
      DF = ZERO
      DDF = ZERO
      DO 30 I = 1, 3
         TEMP = ONE / ( DSCALE( I )-TAU )
         TEMP1 = ZSCALE( I )*TEMP
         TEMP2 = TEMP1*TEMP
         TEMP3 = TEMP2*TEMP
         FC = FC + TEMP1 / DSCALE( I )
         DF = DF + TEMP2
         DDF = DDF + TEMP3
30    CONTINUE
      F = FINIT + TAU*FC
!
      IF( ABS( F ).LE.ZERO ) GO TO 60
      IF( F .LE. ZERO )THEN
         LBD = TAU
      ELSE
         UBD = TAU
      END IF
!
!        Iteration begins -- Use Gragg-Thornton-Warner cubic convergent
!                            scheme
!
!     It is not hard to see that
!
!           1) Iterations will go up monotonically
!              if FINIT < 0;
!
!           2) Iterations will go down monotonically
!              if FINIT > 0.
!
      ITER = NITER + 1
!
      DO 50 NITER = ITER, MAXIT
!
         IF( ORGATI ) THEN
            TEMP1 = DSCALE( 2 ) - TAU
            TEMP2 = DSCALE( 3 ) - TAU
         ELSE
            TEMP1 = DSCALE( 1 ) - TAU
            TEMP2 = DSCALE( 2 ) - TAU
         END IF
         A = ( TEMP1+TEMP2 )*F - TEMP1*TEMP2*DF
         B = TEMP1*TEMP2*F
         C = F - ( TEMP1+TEMP2 )*DF + TEMP1*TEMP2*DDF
         TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
         A = A / TEMP
         B = B / TEMP
         C = C / TEMP
         IF( C.EQ.ZERO ) THEN
            ETA = B / A
         ELSE IF( A.LE.ZERO ) THEN
            ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
         IF( F*ETA.GE.ZERO ) THEN
            ETA = -F / DF
         END IF
!
         TAU = TAU + ETA
         IF( TAU .LT. LBD .OR. TAU .GT. UBD ) TAU = ( LBD + UBD )/TWO 
!
         FC = ZERO
         ERRETM = ZERO
         DF = ZERO
         DDF = ZERO
         DO 40 I = 1, 3
            IF ( ( DSCALE( I )-TAU ).NE.ZERO ) THEN
               TEMP = ONE / ( DSCALE( I )-TAU )
               TEMP1 = ZSCALE( I )*TEMP
               TEMP2 = TEMP1*TEMP
               TEMP3 = TEMP2*TEMP
               TEMP4 = TEMP1 / DSCALE( I )
               FC = FC + TEMP4
               ERRETM = ERRETM + ABS( TEMP4 )
               DF = DF + TEMP2
               DDF = DDF + TEMP3
            ELSE
               GO TO 60
            END IF
40       CONTINUE
         F = FINIT + TAU*FC
         ERRETM = EIGHT*( ABS( FINIT )+ABS( TAU )*ERRETM ) + ABS( TAU )*DF
         IF( ( ABS( F ).LE.FOUR*EPS*ERRETM ) .OR. &
              ( (UBD-LBD).LE.FOUR*EPS*ABS(TAU) )  )&
              GO TO 60
         IF( F .LE. ZERO )THEN
            LBD = TAU
         ELSE
            UBD = TAU
         END IF
50    CONTINUE
      INFO = 1
60    CONTINUE
!
!     Undo scaling
!
      IF( SCALE ) TAU = TAU*SCLINV
      RETURN
!
!     End of DLAED6
!
   END SUBROUTINE DLAED6
!
!=
!
!> \brief \b DLAEDA used by sstedc. Computes the Z vector determining the rank-one modification of the diagonal matrix. Used when the original matrix is dense.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAEDA + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaeda.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaeda.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaeda.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR,
!                          GIVCOL, GIVNUM, Q, QPTR, Z, ZTEMP, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            CURLVL, CURPBM, INFO, N, TLVLS
!       ..
!       .. Array Arguments ..
!       INTEGER            GIVCOL( 2, * ), GIVPTR( * ), PERM( * ),
!      $                   PRMPTR( * ), QPTR( * )
!       DOUBLE PRECISION   GIVNUM( 2, * ), Q( * ), Z( * ), ZTEMP( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAEDA computes the Z vector corresponding to the merge step in the
!> CURLVLth step of the merge process with TLVLS steps for the CURPBMth
!> problem.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] TLVLS
!> \verbatim
!>          TLVLS is INTEGER
!>         The total number of merging levels in the overall divide and
!>         conquer tree.
!> \endverbatim
!>
!> \param[in] CURLVL
!> \verbatim
!>          CURLVL is INTEGER
!>         The current level in the overall merge routine,
!>         0 <= curlvl <= tlvls.
!> \endverbatim
!>
!> \param[in] CURPBM
!> \verbatim
!>          CURPBM is INTEGER
!>         The current problem in the current level in the overall
!>         merge routine (counting from upper left to lower right).
!> \endverbatim
!>
!> \param[in] PRMPTR
!> \verbatim
!>          PRMPTR is INTEGER array, dimension (N lg N)
!>         Contains a list of pointers which indicate where in PERM a
!>         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
!>         indicates the size of the permutation and incidentally the
!>         size of the full, non-deflated problem.
!> \endverbatim
!>
!> \param[in] PERM
!> \verbatim
!>          PERM is INTEGER array, dimension (N lg N)
!>         Contains the permutations (from deflation and sorting) to be
!>         applied to each eigenblock.
!> \endverbatim
!>
!> \param[in] GIVPTR
!> \verbatim
!>          GIVPTR is INTEGER array, dimension (N lg N)
!>         Contains a list of pointers which indicate where in GIVCOL a
!>         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
!>         indicates the number of Givens rotations.
!> \endverbatim
!>
!> \param[in] GIVCOL
!> \verbatim
!>          GIVCOL is INTEGER array, dimension (2, N lg N)
!>         Each pair of numbers indicates a pair of columns to take place
!>         in a Givens rotation.
!> \endverbatim
!>
!> \param[in] GIVNUM
!> \verbatim
!>          GIVNUM is DOUBLE PRECISION array, dimension (2, N lg N)
!>         Each number indicates the S value to be used in the
!>         corresponding Givens rotation.
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (N**2)
!>         Contains the square eigenblocks from previous levels, the
!>         starting positions for blocks are given by QPTR.
!> \endverbatim
!>
!> \param[in] QPTR
!> \verbatim
!>          QPTR is INTEGER array, dimension (N+2)
!>         Contains a list of pointers which indicate where in Q an
!>         eigenblock is stored.  SQRT( QPTR(i+1) - QPTR(i) ) indicates
!>         the size of the block.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (N)
!>         On output this vector contains the updating vector (the last
!>         row of the first sub-eigenvector matrix and the first row of
!>         the second sub-eigenvector matrix).
!> \endverbatim
!>
!> \param[out] ZTEMP
!> \verbatim
!>          ZTEMP is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
   SUBROUTINE DLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR,&
        GIVCOL, GIVNUM, Q, QPTR, Z, ZTEMP, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
     INTEGER            CURLVL, CURPBM, INFO, N, TLVLS
!     ..
!     .. Array Arguments ..
     INTEGER            GIVCOL( 2, * ), GIVPTR( * ), PERM( * ),&
          PRMPTR( * ), QPTR( * )
     DOUBLE PRECISION   GIVNUM( 2, * ), Q( * ), Z( * ), ZTEMP( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
     DOUBLE PRECISION   ZERO, HALF, ONE
     PARAMETER          ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
     INTEGER            BSIZ1, BSIZ2, CURR, I, K, MID, PSIZ1, PSIZ2,&
          PTR, ZPTR1
!     ..
!     .. External Subroutines ..
!     EXTERNAL           DCOPY, DGEMV, DROT, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          DBLE, INT, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
     INFO = 0
!
     IF( N.LT.0 ) THEN
        INFO = -1
     END IF
     IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLAEDA', -INFO )
        RETURN
     END IF
!
!     Quick return if possible
!
     IF( N.EQ.0 ) RETURN
!
!     Determine location of first number in second half.
!
     MID = N / 2 + 1
!
!     Gather last/first rows of appropriate eigenblocks into center of Z
!
     PTR = 1
!
!     Determine location of lowest level subproblem in the full storage
!     scheme
!
     CURR = PTR + CURPBM*2**CURLVL + 2**( CURLVL-1 ) - 1
!
!     Determine size of these matrices.  We add HALF to the value of
!     the SQRT in case the machine underestimates one of these square
!     roots.
!
     BSIZ1 = INT( HALF+SQRT( DBLE( QPTR( CURR+1 )-QPTR( CURR ) ) ) )
     BSIZ2 = INT( HALF+SQRT( DBLE( QPTR( CURR+2 )-QPTR( CURR+1 ) ) ) )
     DO 10 K = 1, MID - BSIZ1 - 1
        Z( K ) = ZERO
10   CONTINUE
     CALL DCOPY( BSIZ1, Q( QPTR( CURR )+BSIZ1-1 ), BSIZ1,&
          Z( MID-BSIZ1 ), 1 )
     CALL DCOPY( BSIZ2, Q( QPTR( CURR+1 ) ), BSIZ2, Z( MID ), 1 )
     DO 20 K = MID + BSIZ2, N
        Z( K ) = ZERO
20   CONTINUE
!
!     Loop through remaining levels 1 -> CURLVL applying the Givens
!     rotations and permutation and then multiplying the center matrices
!     against the current Z.
!
     PTR = 2**TLVLS + 1
     DO 70 K = 1, CURLVL - 1
        CURR = PTR + CURPBM*2**( CURLVL-K ) + 2**( CURLVL-K-1 ) - 1
        PSIZ1 = PRMPTR( CURR+1 ) - PRMPTR( CURR )
        PSIZ2 = PRMPTR( CURR+2 ) - PRMPTR( CURR+1 )
        ZPTR1 = MID - PSIZ1
!
!       Apply Givens at CURR and CURR+1
!
        DO 30 I = GIVPTR( CURR ), GIVPTR( CURR+1 ) - 1
           CALL DROT( 1, Z( ZPTR1+GIVCOL( 1, I )-1 ), 1,&
                Z( ZPTR1+GIVCOL( 2, I )-1 ), 1, GIVNUM( 1, I ),&
                GIVNUM( 2, I ) )
30      CONTINUE
        DO 40 I = GIVPTR( CURR+1 ), GIVPTR( CURR+2 ) - 1
           CALL DROT( 1, Z( MID-1+GIVCOL( 1, I ) ), 1,&
                Z( MID-1+GIVCOL( 2, I ) ), 1, GIVNUM( 1, I ),&
                GIVNUM( 2, I ) )
40      CONTINUE
        PSIZ1 = PRMPTR( CURR+1 ) - PRMPTR( CURR )
        PSIZ2 = PRMPTR( CURR+2 ) - PRMPTR( CURR+1 )
        DO 50 I = 0, PSIZ1 - 1
           ZTEMP( I+1 ) = Z( ZPTR1+PERM( PRMPTR( CURR )+I )-1 )
50      CONTINUE
        DO 60 I = 0, PSIZ2 - 1
           ZTEMP( PSIZ1+I+1 ) = Z( MID+PERM( PRMPTR( CURR+1 )+I )-1 )
60      CONTINUE
!
!        Multiply Blocks at CURR and CURR+1
!
!        Determine size of these matrices.  We add HALF to the value of
!        the SQRT in case the machine underestimates one of these
!        square roots.
!
        BSIZ1 = INT( HALF+SQRT( DBLE( QPTR( CURR+1 )-QPTR( CURR ) ) ) )
        BSIZ2 = INT( HALF+SQRT( DBLE( QPTR( CURR+2 )-QPTR( CURR+&
             1 ) ) ) )
        IF( BSIZ1.GT.0 ) THEN
           CALL DGEMV( 'T', BSIZ1, BSIZ1, ONE, Q( QPTR( CURR ) ),&
                BSIZ1, ZTEMP( 1 ), 1, ZERO, Z( ZPTR1 ), 1 )
        END IF
        CALL DCOPY( PSIZ1-BSIZ1, ZTEMP( BSIZ1+1 ), 1, Z( ZPTR1+BSIZ1 ), 1 )
        IF( BSIZ2.GT.0 ) THEN
           CALL DGEMV( 'T', BSIZ2, BSIZ2, ONE, Q( QPTR( CURR+1 ) ),&
                BSIZ2, ZTEMP( PSIZ1+1 ), 1, ZERO, Z( MID ), 1 )
        END IF
        CALL DCOPY( PSIZ2-BSIZ2, ZTEMP( PSIZ1+BSIZ2+1 ), 1,&
             Z( MID+BSIZ2 ), 1 )
!
        PTR = PTR + 2**( TLVLS-K )
70   CONTINUE
!
     RETURN
!
!     End of DLAEDA
!
  END SUBROUTINE DLAEDA
!
!=
!
!> \brief \b DOPGTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DOPGTR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dopgtr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dopgtr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dopgtr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DOPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDQ, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AP( * ), Q( LDQ, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DOPGTR generates a real orthogonal matrix Q which is defined as the
!> product of n-1 elementary reflectors H(i) of order n, as returned by
!> DSPTRD using packed storage:
!>
!> if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!>
!> if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U': Upper triangular packed storage used in previous
!>                 call to DSPTRD;
!>          = 'L': Lower triangular packed storage used in previous
!>                 call to DSPTRD.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix Q. N >= 0.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          The vectors which define the elementary reflectors, as
!>          returned by DSPTRD.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (N-1)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DSPTRD.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!>          The N-by-N orthogonal matrix Q.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q. LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N-1)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
  SUBROUTINE DOPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    CHARACTER          UPLO
    INTEGER            INFO, LDQ, N
!     ..
!     .. Array Arguments ..
    DOUBLE PRECISION   AP( * ), Q( LDQ, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
    DOUBLE PRECISION   ZERO, ONE
    PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
    LOGICAL            UPPER
    INTEGER            I, IINFO, IJ, J
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DORG2L, DORG2R, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
    INFO = 0
    UPPER = LSAME( UPLO, 'U' )
    IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
       INFO = -1
    ELSE IF( N.LT.0 ) THEN
       INFO = -2
    ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
       INFO = -6
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DOPGTR', -INFO )
       RETURN
    END IF
!
!     Quick return if possible
!
    IF( N.EQ.0 ) RETURN
!
    IF( UPPER ) THEN
!
!        Q was determined by a call to DSPTRD with UPLO = 'U'
!
!        Unpack the vectors which define the elementary reflectors and
!        set the last row and column of Q equal to those of the unit
!        matrix
!
       IJ = 2
       DO 20 J = 1, N - 1
          DO 10 I = 1, J - 1
             Q( I, J ) = AP( IJ )
             IJ = IJ + 1
10        CONTINUE
          IJ = IJ + 2
          Q( N, J ) = ZERO
20     CONTINUE
       DO 30 I = 1, N - 1
          Q( I, N ) = ZERO
30     CONTINUE
       Q( N, N ) = ONE
!
!        Generate Q(1:n-1,1:n-1)
!
       CALL DORG2L( N-1, N-1, N-1, Q, LDQ, TAU, WORK, IINFO )
!
    ELSE
!
!        Q was determined by a call to DSPTRD with UPLO = 'L'.
!
!        Unpack the vectors which define the elementary reflectors and
!        set the first row and column of Q equal to those of the unit
!        matrix
!
       Q( 1, 1 ) = ONE
       DO 40 I = 2, N
          Q( I, 1 ) = ZERO
40     CONTINUE
       IJ = 3
       DO 60 J = 2, N
          Q( 1, J ) = ZERO
          DO 50 I = J + 1, N
             Q( I, J ) = AP( IJ )
             IJ = IJ + 1
50        CONTINUE
          IJ = IJ + 2
60     CONTINUE
       IF( N.GT.1 ) THEN
!
!           Generate Q(2:n,2:n)
!
          CALL DORG2R( N-1, N-1, N-1, Q( 2, 2 ), LDQ, TAU, WORK, IINFO )
       END IF
    END IF
    RETURN
!
!     End of DOPGTR
!
  END SUBROUTINE DOPGTR
!
!=
!
!> \brief \b DORG2L generates all or part of the orthogonal matrix Q from a QL factorization determined by sgeqlf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DORG2L + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2l.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2l.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2l.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORG2L generates an m by n real matrix Q with orthonormal columns,
!> which is defined as the last n columns of a product of k elementary
!> reflectors of order m
!>
!>       Q  =  H(k) . . . H(2) H(1)
!>
!> as returned by DGEQLF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the (n-k+i)-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by DGEQLF in the last k columns of its array
!>          argument A.
!>          On exit, the m by n matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGEQLF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument has an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
  SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
    INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
    DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
    DOUBLE PRECISION   ONE, ZERO
    PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
    INTEGER            I, II, J, L
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLARF, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
    INFO = 0
    IF( M.LT.0 ) THEN
       INFO = -1
    ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
       INFO = -2
    ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
       INFO = -3
    ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
       INFO = -5
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DORG2L', -INFO )
       RETURN
    END IF
!
!     Quick return if possible
!
    IF( N.LE.0 ) RETURN
!
!     Initialise columns 1:n-k to columns of the unit matrix
!
    DO 20 J = 1, N - K
       DO 10 L = 1, M
          A( L, J ) = ZERO
10     CONTINUE
       A( M-N+J, J ) = ONE
20  CONTINUE
!
    DO 40 I = 1, K
       II = N - K + I
!
!        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
!
       A( M-N+II, II ) = ONE
       CALL DLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A,&
            LDA, WORK )
       CALL DSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
       A( M-N+II, II ) = ONE - TAU( I )
!
!        Set A(m-k+i+1:m,n-k+i) to zero
!
       DO 30 L = M - N + II + 1, M
          A( L, II ) = ZERO
30     CONTINUE
40  CONTINUE
    RETURN
!
!     End of DORG2L
!
  END SUBROUTINE DORG2L
!
!=
!
!> \brief \b DORG2R generates all or part of the orthogonal matrix Q from a QR factorization determined by sgeqrf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DORG2R + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2r.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2r.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2r.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORG2R generates an m by n real matrix Q with orthonormal columns,
!> which is defined as the first n columns of a product of k elementary
!> reflectors of order m
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by DGEQRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the i-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by DGEQRF in the first k columns of its array
!>          argument A.
!>          On exit, the m-by-n matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGEQRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument has an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
  SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
    INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
    DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
    DOUBLE PRECISION   ONE, ZERO
    PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
    INTEGER            I, J, L
!     ..
!     .. External Subroutines ..
!    EXTERNAL           DLARF, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
    INFO = 0
    IF( M.LT.0 ) THEN
       INFO = -1
    ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
       INFO = -2
    ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
       INFO = -3
    ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
       INFO = -5
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DORG2R', -INFO )
       RETURN
    END IF
!
!     Quick return if possible
!
    IF( N.LE.0 ) RETURN
!
!     Initialise columns k+1:n to columns of the unit matrix
!
    DO 20 J = K + 1, N
       DO 10 L = 1, M
          A( L, J ) = ZERO
10     CONTINUE
       A( J, J ) = ONE
20  CONTINUE
!
    DO 40 I = K, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the left
!
       IF( I.LT.N ) THEN
          A( I, I ) = ONE
          CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),&
               A( I, I+1 ), LDA, WORK )
       END IF
       IF( I.LT.M ) CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
       A( I, I ) = ONE - TAU( I )
!
!        Set A(1:i-1,i) to zero
!
       DO 30 L = 1, I - 1
          A( L, I ) = ZERO
30     CONTINUE
40  CONTINUE
    RETURN
!
!     End of DORG2R
!
  END SUBROUTINE DORG2R
!
!========================================================================
!
! BLAS
!
!> \brief \b DROT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION C,S
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DROT applies a plane rotation.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
  SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    DOUBLE PRECISION C,S
    INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
    DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
    DOUBLE PRECISION DTEMP
    INTEGER I,IX,IY
!     ..
    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!       code for both increments equal to 1
!
       DO I = 1,N
          DTEMP = C*DX(I) + S*DY(I)
          DY(I) = C*DY(I) - S*DX(I)
          DX(I) = DTEMP
       END DO
    ELSE
!
!       code for unequal increments or equal increments not equal
!         to 1
!
       IX = 1
       IY = 1
       IF (INCX.LT.0) IX = (-N+1)*INCX + 1
       IF (INCY.LT.0) IY = (-N+1)*INCY + 1
       DO I = 1,N
          DTEMP = C*DX(IX) + S*DY(IY)
          DY(IY) = C*DY(IY) - S*DX(IX)
          DX(IX) = DTEMP
          IX = IX + INCX
          IY = IY + INCY
       END DO
    END IF
    RETURN
  END SUBROUTINE DROT
!
!=
!
!> \brief \b DNRM2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
! 
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION X(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DNRM2 returns the euclidean norm of a vector via the function
!> name, so that
!>
!>    DNRM2 := sqrt( x'*x )
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  -- This version written on 25-October-1982.
!>     Modified on 14-October-1993 to inline the call to DLASSQ.
!>     Sven Hammarling, Nag Ltd.
!> \endverbatim
!>
!  =====================================================================
  DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    INTEGER INCX,N
!     ..
!     .. Array Arguments ..
    DOUBLE PRECISION X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
    DOUBLE PRECISION ONE,ZERO
    PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
    DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
    INTEGER IX
!     ..
!     .. Intrinsic Functions ..
    INTRINSIC ABS,SQRT
!     ..
    IF (N.LT.1 .OR. INCX.LT.1) THEN
       NORM = ZERO
    ELSE IF (N.EQ.1) THEN
       NORM = ABS(X(1))
    ELSE
       SCALE = ZERO
       SSQ = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
!
       DO 10 IX = 1,1 + (N-1)*INCX,INCX
          IF (X(IX).NE.ZERO) THEN
             ABSXI = ABS(X(IX))
             IF (SCALE.LT.ABSXI) THEN
                SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                SCALE = ABSXI
             ELSE
                SSQ = SSQ + (ABSXI/SCALE)**2
             END IF
          END IF
10     CONTINUE
       NORM = SCALE*SQRT(SSQ)
    END IF
!
    DNRM2 = NORM
    RETURN
!
!     End of DNRM2.
!
  END FUNCTION DNRM2
!
!=
!
!> \brief \b DGER
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA
!       INTEGER INCX,INCY,LDA,M,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGER   performs the rank 1 operation
!>
!>    A := alpha*x*y**T + A,
!>
!> where alpha is a scalar, x is an m element vector, y is an n element
!> vector and A is an m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of the matrix A.
!>           M must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array of dimension at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the m
!>           element vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array of dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients. On exit, A is
!>           overwritten by the updated matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, m ).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
  SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!  -- Reference BLAS level2 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    DOUBLE PRECISION ALPHA
    INTEGER INCX,INCY,LDA,M,N
!     ..
!     .. Array Arguments ..
    DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
    DOUBLE PRECISION ZERO
    PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
    DOUBLE PRECISION TEMP
    INTEGER I,INFO,IX,J,JY,KX
!     ..
!     .. External Subroutines ..
!    EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!    INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
    INFO = 0
    IF (M.LT.0) THEN
       INFO = 1
    ELSE IF (N.LT.0) THEN
       INFO = 2
    ELSE IF (INCX.EQ.0) THEN
       INFO = 5
    ELSE IF (INCY.EQ.0) THEN
       INFO = 7
    ELSE IF (LDA.LT.MAX(1,M)) THEN
       INFO = 9
    END IF
    IF (INFO.NE.0) THEN
       CALL XERBLA('DGER  ',INFO)
       RETURN
    END IF
!
!     Quick return if possible.
!
    IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
    IF (INCY.GT.0) THEN
       JY = 1
    ELSE
       JY = 1 - (N-1)*INCY
    END IF
    IF (INCX.EQ.1) THEN
       DO 20 J = 1,N
          IF (Y(JY).NE.ZERO) THEN
             TEMP = ALPHA*Y(JY)
             DO 10 I = 1,M
                A(I,J) = A(I,J) + X(I)*TEMP
10           CONTINUE
          END IF
          JY = JY + INCY
20     CONTINUE
    ELSE
       IF (INCX.GT.0) THEN
          KX = 1
       ELSE
          KX = 1 - (M-1)*INCX
       END IF
       DO 40 J = 1,N
          IF (Y(JY).NE.ZERO) THEN
             TEMP = ALPHA*Y(JY)
             IX = KX
             DO 30 I = 1,M
                A(I,J) = A(I,J) + X(IX)*TEMP
                IX = IX + INCX
30           CONTINUE
          END IF
          JY = JY + INCY
40     CONTINUE
      END IF
!
   RETURN
!
!     End of DGER  .
!
 END SUBROUTINE DGER
!
!=
!
!> \brief \b DSPMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER INCX,INCY,N
!       CHARACTER UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION AP(*),X(*),Y(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSPMV  performs the matrix-vector operation
!>
!>    y := alpha*A*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are n element vectors and
!> A is an n by n symmetric matrix, supplied in packed form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the matrix A is supplied in the packed
!>           array AP as follows:
!>
!>              UPLO = 'U' or 'u'   The upper triangular part of A is
!>                                  supplied in AP.
!>
!>              UPLO = 'L' or 'l'   The lower triangular part of A is
!>                                  supplied in AP.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array of DIMENSION at least
!>           ( ( n*( n + 1 ) )/2 ).
!>           Before entry with UPLO = 'U' or 'u', the array AP must
!>           contain the upper triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!>           and a( 2, 2 ) respectively, and so on.
!>           Before entry with UPLO = 'L' or 'l', the array AP must
!>           contain the lower triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!>           and a( 3, 1 ) respectively, and so on.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array of dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
!>           element vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION.
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array of dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y. On exit, Y is overwritten by the updated
!>           vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
 SUBROUTINE DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
!
!  -- Reference BLAS level2 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
   DOUBLE PRECISION ALPHA,BETA
   INTEGER INCX,INCY,N
   CHARACTER UPLO
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION AP(*),X(*),Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION ONE,ZERO
   PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
   DOUBLE PRECISION TEMP1,TEMP2
   INTEGER I,INFO,IX,IY,J,JX,JY,K,KK,KX,KY
!     ..
!     .. External Functions ..
!   LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!
!     Test the input parameters.
!
   INFO = 0
   IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
      INFO = 1
   ELSE IF (N.LT.0) THEN
      INFO = 2
   ELSE IF (INCX.EQ.0) THEN
      INFO = 6
   ELSE IF (INCY.EQ.0) THEN
      INFO = 9
   END IF
   IF (INFO.NE.0) THEN
      CALL XERBLA('DSPMV ',INFO)
      RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!
!     Set up the start points in  X  and  Y.
!
   IF (INCX.GT.0) THEN
      KX = 1
   ELSE
      KX = 1 - (N-1)*INCX
   END IF
   IF (INCY.GT.0) THEN
      KY = 1
   ELSE
      KY = 1 - (N-1)*INCY
   END IF
!
!     Start the operations. In this version the elements of the array AP
!     are accessed sequentially with one pass through AP.
!
!     First form  y := beta*y.
!
   IF (BETA.NE.ONE) THEN
      IF (INCY.EQ.1) THEN
         IF (BETA.EQ.ZERO) THEN
            DO 10 I = 1,N
               Y(I) = ZERO
10          CONTINUE
         ELSE
            DO 20 I = 1,N
               Y(I) = BETA*Y(I)
20          CONTINUE
         END IF
      ELSE
         IY = KY
         IF (BETA.EQ.ZERO) THEN
            DO 30 I = 1,N
               Y(IY) = ZERO
               IY = IY + INCY
30          CONTINUE
         ELSE
            DO 40 I = 1,N
               Y(IY) = BETA*Y(IY)
               IY = IY + INCY
40          CONTINUE
         END IF
      END IF
   END IF
   IF (ALPHA.EQ.ZERO) RETURN
   KK = 1
   IF (LSAME(UPLO,'U')) THEN
!
!        Form  y  when AP contains the upper triangle.
!
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
         DO 60 J = 1,N
            TEMP1 = ALPHA*X(J)
            TEMP2 = ZERO
            K = KK
            DO 50 I = 1,J - 1
               Y(I) = Y(I) + TEMP1*AP(K)
               TEMP2 = TEMP2 + AP(K)*X(I)
               K = K + 1
50          CONTINUE
            Y(J) = Y(J) + TEMP1*AP(KK+J-1) + ALPHA*TEMP2
            KK = KK + J
60       CONTINUE
      ELSE
         JX = KX
         JY = KY
         DO 80 J = 1,N
            TEMP1 = ALPHA*X(JX)
            TEMP2 = ZERO
            IX = KX
            IY = KY
            DO 70 K = KK,KK + J - 2
               Y(IY) = Y(IY) + TEMP1*AP(K)
               TEMP2 = TEMP2 + AP(K)*X(IX)
               IX = IX + INCX
               IY = IY + INCY
70          CONTINUE
            Y(JY) = Y(JY) + TEMP1*AP(KK+J-1) + ALPHA*TEMP2
            JX = JX + INCX
            JY = JY + INCY
            KK = KK + J
80       CONTINUE
      END IF
   ELSE
!
!        Form  y  when AP contains the lower triangle.
!
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
         DO 100 J = 1,N
            TEMP1 = ALPHA*X(J)
            TEMP2 = ZERO
            Y(J) = Y(J) + TEMP1*AP(KK)
            K = KK + 1
            DO 90 I = J + 1,N
               Y(I) = Y(I) + TEMP1*AP(K)
               TEMP2 = TEMP2 + AP(K)*X(I)
               K = K + 1
90          CONTINUE
            Y(J) = Y(J) + ALPHA*TEMP2
            KK = KK + (N-J+1)
100      CONTINUE
      ELSE
         JX = KX
         JY = KY
         DO 120 J = 1,N
            TEMP1 = ALPHA*X(JX)
            TEMP2 = ZERO
            Y(JY) = Y(JY) + TEMP1*AP(KK)
            IX = JX
            IY = JY
            DO 110 K = KK + 1,KK + N - J
               IX = IX + INCX
               IY = IY + INCY
               Y(IY) = Y(IY) + TEMP1*AP(K)
               TEMP2 = TEMP2 + AP(K)*X(IX)
110         CONTINUE
            Y(JY) = Y(JY) + ALPHA*TEMP2
            JX = JX + INCX
            JY = JY + INCY
            KK = KK + (N-J+1)
120      CONTINUE
      END IF
   END IF
!
   RETURN
!
!     End of DSPMV .
!
 END SUBROUTINE DSPMV
!
!=
!
!> \brief \b DAXPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION DA
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DAXPY constant times a vector plus a vector.
!>    uses unrolled loops for increments equal to one.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
 SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
   DOUBLE PRECISION DA
   INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
   INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC MOD
!     ..
   IF (N.LE.0) RETURN
   IF (DA.EQ.0.0d0) RETURN
   IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
      M = MOD(N,4)
      IF (M.NE.0) THEN
         DO I = 1,M
            DY(I) = DY(I) + DA*DX(I)
         END DO
      END IF
      IF (N.LT.4) RETURN
      MP1 = M + 1
      DO I = MP1,N,4
         DY(I) = DY(I) + DA*DX(I)
         DY(I+1) = DY(I+1) + DA*DX(I+1)
         DY(I+2) = DY(I+2) + DA*DX(I+2)
         DY(I+3) = DY(I+3) + DA*DX(I+3)
      END DO
   ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO I = 1,N
         DY(IY) = DY(IY) + DA*DX(IX)
         IX = IX + INCX
         IY = IY + INCY
      END DO
   END IF
   RETURN
 END SUBROUTINE DAXPY
!
!=
!
!> \brief \b DSPR2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA
!       INTEGER INCX,INCY,N
!       CHARACTER UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION AP(*),X(*),Y(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSPR2  performs the symmetric rank 2 operation
!>
!>    A := alpha*x*y**T + alpha*y*x**T + A,
!>
!> where alpha is a scalar, x and y are n element vectors and A is an
!> n by n symmetric matrix, supplied in packed form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the matrix A is supplied in the packed
!>           array AP as follows:
!>
!>              UPLO = 'U' or 'u'   The upper triangular part of A is
!>                                  supplied in AP.
!>
!>              UPLO = 'L' or 'l'   The lower triangular part of A is
!>                                  supplied in AP.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array of dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
!>           element vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array of dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!>
!> \param[in,out] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array of DIMENSION at least
!>           ( ( n*( n + 1 ) )/2 ).
!>           Before entry with  UPLO = 'U' or 'u', the array AP must
!>           contain the upper triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!>           and a( 2, 2 ) respectively, and so on. On exit, the array
!>           AP is overwritten by the upper triangular part of the
!>           updated matrix.
!>           Before entry with UPLO = 'L' or 'l', the array AP must
!>           contain the lower triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!>           and a( 3, 1 ) respectively, and so on. On exit, the array
!>           AP is overwritten by the lower triangular part of the
!>           updated matrix.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
 SUBROUTINE DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
!
!  -- Reference BLAS level2 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
   DOUBLE PRECISION ALPHA
   INTEGER INCX,INCY,N
   CHARACTER UPLO
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION AP(*),X(*),Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION ZERO
   PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
   DOUBLE PRECISION TEMP1,TEMP2
   INTEGER I,INFO,IX,IY,J,JX,JY,K,KK,KX,KY
!     ..
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!
!     Test the input parameters.
!
   INFO = 0
   IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
      INFO = 1
   ELSE IF (N.LT.0) THEN
      INFO = 2
   ELSE IF (INCX.EQ.0) THEN
      INFO = 5
   ELSE IF (INCY.EQ.0) THEN
      INFO = 7
   END IF
   IF (INFO.NE.0) THEN
      CALL XERBLA('DSPR2 ',INFO)
      RETURN
   END IF
!
!     Quick return if possible.
!
   IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!
!     Set up the start points in X and Y if the increments are not both
!     unity.
!
   IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
      IF (INCX.GT.0) THEN
         KX = 1
      ELSE
         KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
         KY = 1
      ELSE
         KY = 1 - (N-1)*INCY
      END IF
      JX = KX
      JY = KY
   END IF
!
!     Start the operations. In this version the elements of the array AP
!     are accessed sequentially with one pass through AP.
!
   KK = 1
   IF (LSAME(UPLO,'U')) THEN
!
!        Form  A  when upper triangle is stored in AP.
!
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
         DO 20 J = 1,N
            IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
               TEMP1 = ALPHA*Y(J)
               TEMP2 = ALPHA*X(J)
               K = KK
               DO 10 I = 1,J
                  AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                  K = K + 1
10             CONTINUE
            END IF
            KK = KK + J
20       CONTINUE
      ELSE
         DO 40 J = 1,N
            IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
               TEMP1 = ALPHA*Y(JY)
               TEMP2 = ALPHA*X(JX)
               IX = KX
               IY = KY
               DO 30 K = KK,KK + J - 1
                  AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                  IX = IX + INCX
                  IY = IY + INCY
30             CONTINUE
            END IF
            JX = JX + INCX
            JY = JY + INCY
            KK = KK + J
40       CONTINUE
      END IF
   ELSE
!
!        Form  A  when lower triangle is stored in AP.
!
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
         DO 60 J = 1,N
            IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
               TEMP1 = ALPHA*Y(J)
               TEMP2 = ALPHA*X(J)
               K = KK
               DO 50 I = J,N
                  AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                  K = K + 1
50             CONTINUE
            END IF
            KK = KK + N - J + 1
60       CONTINUE
      ELSE
         DO 80 J = 1,N
            IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
               TEMP1 = ALPHA*Y(JY)
               TEMP2 = ALPHA*X(JX)
               IX = JX
               IY = JY
               DO 70 K = KK,KK + N - J
                  AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                  IX = IX + INCX
                  IY = IY + INCY
70             CONTINUE
            END IF
            JX = JX + INCX
            JY = JY + INCY
            KK = KK + N - J + 1
80       CONTINUE
      END IF
   END IF
!
   RETURN
!
!     End of DSPR2 .
!
 END SUBROUTINE DSPR2
!
!=
!
!> \brief \b DTRMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
! 
!       .. Scalar Arguments ..
!       INTEGER INCX,LDA,N
!       CHARACTER DIAG,TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),X(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTRMV  performs one of the matrix-vector operations
!>
!>    x := A*x,   or   x := A**T*x,
!>
!> where x is an n element vector and  A is an n by n unit, or non-unit,
!> upper or lower triangular matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the matrix is an upper or
!>           lower triangular matrix as follows:
!>
!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!>
!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   x := A*x.
!>
!>              TRANS = 'T' or 't'   x := A**T*x.
!>
!>              TRANS = 'C' or 'c'   x := A**T*x.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not A is unit
!>           triangular as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!>           Before entry with  UPLO = 'U' or 'u', the leading n by n
!>           upper triangular part of the array A must contain the upper
!>           triangular matrix and the strictly lower triangular part of
!>           A is not referenced.
!>           Before entry with UPLO = 'L' or 'l', the leading n by n
!>           lower triangular part of the array A must contain the lower
!>           triangular matrix and the strictly upper triangular part of
!>           A is not referenced.
!>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!>           A are not referenced either, but are assumed to be unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array of dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
!>           element vector x. On exit, X is overwritten with the
!>           tranformed vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
 SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!
!  -- Reference BLAS level2 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
   INTEGER INCX,LDA,N
   CHARACTER DIAG,TRANS,UPLO
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION A(LDA,*),X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
   DOUBLE PRECISION ZERO
   PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
   DOUBLE PRECISION TEMP
   INTEGER I,INFO,IX,J,JX,KX
   LOGICAL NOUNIT
!     ..
!     .. External Functions ..
!   LOGICAL LSAME
!   EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
   INFO = 0
   IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
      INFO = 1
   ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.&
        .NOT.LSAME(TRANS,'C')) THEN
      INFO = 2
   ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
      INFO = 3
   ELSE IF (N.LT.0) THEN
      INFO = 4
   ELSE IF (LDA.LT.MAX(1,N)) THEN
      INFO = 6
   ELSE IF (INCX.EQ.0) THEN
      INFO = 8
   END IF
   IF (INFO.NE.0) THEN
      CALL XERBLA('DTRMV ',INFO)
      RETURN
   END IF
!
!     Quick return if possible.
!
   IF (N.EQ.0) RETURN
!
   NOUNIT = LSAME(DIAG,'N')
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
   IF (INCX.LE.0) THEN
      KX = 1 - (N-1)*INCX
   ELSE IF (INCX.NE.1) THEN
      KX = 1
   END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
   IF (LSAME(TRANS,'N')) THEN
!
!        Form  x := A*x.
!
      IF (LSAME(UPLO,'U')) THEN
         IF (INCX.EQ.1) THEN
            DO 20 J = 1,N
               IF (X(J).NE.ZERO) THEN
                  TEMP = X(J)
                  DO 10 I = 1,J - 1
                     X(I) = X(I) + TEMP*A(I,J)
10                CONTINUE
                  IF (NOUNIT) X(J) = X(J)*A(J,J)
               END IF
20          CONTINUE
         ELSE
            JX = KX
            DO 40 J = 1,N
               IF (X(JX).NE.ZERO) THEN
                  TEMP = X(JX)
                  IX = KX
                  DO 30 I = 1,J - 1
                     X(IX) = X(IX) + TEMP*A(I,J)
                     IX = IX + INCX
30                CONTINUE
                  IF (NOUNIT) X(JX) = X(JX)*A(J,J)
               END IF
               JX = JX + INCX
40          CONTINUE
         END IF
      ELSE
         IF (INCX.EQ.1) THEN
            DO 60 J = N,1,-1
               IF (X(J).NE.ZERO) THEN
                  TEMP = X(J)
                  DO 50 I = N,J + 1,-1
                     X(I) = X(I) + TEMP*A(I,J)
50                CONTINUE
                  IF (NOUNIT) X(J) = X(J)*A(J,J)
               END IF
60          CONTINUE
         ELSE
            KX = KX + (N-1)*INCX
            JX = KX
            DO 80 J = N,1,-1
               IF (X(JX).NE.ZERO) THEN
                  TEMP = X(JX)
                  IX = KX
                  DO 70 I = N,J + 1,-1
                     X(IX) = X(IX) + TEMP*A(I,J)
                     IX = IX - INCX
70                CONTINUE
                  IF (NOUNIT) X(JX) = X(JX)*A(J,J)
               END IF
               JX = JX - INCX
80          CONTINUE
         END IF
      END IF
   ELSE
!
!        Form  x := A**T*x.
!
      IF (LSAME(UPLO,'U')) THEN
         IF (INCX.EQ.1) THEN
            DO 100 J = N,1,-1
               TEMP = X(J)
               IF (NOUNIT) TEMP = TEMP*A(J,J)
               DO 90 I = J - 1,1,-1
                  TEMP = TEMP + A(I,J)*X(I)
90             CONTINUE
               X(J) = TEMP
100         CONTINUE
         ELSE
            JX = KX + (N-1)*INCX
            DO 120 J = N,1,-1
               TEMP = X(JX)
               IX = JX
               IF (NOUNIT) TEMP = TEMP*A(J,J)
               DO 110 I = J - 1,1,-1
                  IX = IX - INCX
                  TEMP = TEMP + A(I,J)*X(IX)
110            CONTINUE
               X(JX) = TEMP
               JX = JX - INCX
120         CONTINUE
         END IF
      ELSE
         IF (INCX.EQ.1) THEN
            DO 140 J = 1,N
               TEMP = X(J)
               IF (NOUNIT) TEMP = TEMP*A(J,J)
               DO 130 I = J + 1,N
                  TEMP = TEMP + A(I,J)*X(I)
130            CONTINUE
               X(J) = TEMP
140         CONTINUE
         ELSE
            JX = KX
            DO 160 J = 1,N
               TEMP = X(JX)
               IX = JX
               IF (NOUNIT) TEMP = TEMP*A(J,J)
               DO 150 I = J + 1,N
                  IX = IX + INCX
                  TEMP = TEMP + A(I,J)*X(IX)
150            CONTINUE
               X(JX) = TEMP
               JX = JX + INCX
160         CONTINUE
         END IF
      END IF
   END IF
!
   RETURN
!
!     End of DTRMV .
!
 END SUBROUTINE DTRMV
!> \brief \b DTRMM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA
!       INTEGER LDA,LDB,M,N
!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),B(LDB,*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTRMM  performs one of the matrix-matrix operations
!>
!>    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!>
!> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!>
!>    op( A ) = A   or   op( A ) = A**T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>           On entry,  SIDE specifies whether  op( A ) multiplies B from
!>           the left or right as follows:
!>
!>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!>
!>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the matrix A is an upper or
!>           lower triangular matrix as follows:
!>
!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!>
!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>           On entry, TRANSA specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA = 'N' or 'n'   op( A ) = A.
!>
!>              TRANSA = 'T' or 't'   op( A ) = A**T.
!>
!>              TRANSA = 'C' or 'c'   op( A ) = A**T.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not A is unit triangular
!>           as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of B. M must be at
!>           least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of B.  N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!>           zero then  A is not referenced and  B need not be set before
!>           entry.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>           A is DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!>           upper triangular part of the array  A must contain the upper
!>           triangular matrix  and the strictly lower triangular part of
!>           A is not referenced.
!>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!>           lower triangular part of the array  A must contain the lower
!>           triangular matrix  and the strictly upper triangular part of
!>           A is not referenced.
!>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!>           A  are not referenced either,  but are assumed to be  unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!>           then LDA must be at least max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!>           Before entry,  the leading  m by n part of the array  B must
!>           contain the matrix  B,  and  on exit  is overwritten  by the
!>           transformed matrix.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
!>           max( 1, m ).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup double_blas_level3
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 3 Blas routine.
!>
!>  -- Written on 8-February-1989.
!>     Jack Dongarra, Argonne National Laboratory.
!>     Iain Duff, AERE Harwell.
!>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!>     Sven Hammarling, Numerical Algorithms Group Ltd.
!> \endverbatim
!>
!  =====================================================================
 SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!  -- Reference BLAS level3 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
   DOUBLE PRECISION ALPHA
   INTEGER LDA,LDB,M,N
   CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
   DOUBLE PRECISION A(LDA,*),B(LDB,*)
!     ..
!
!  =====================================================================
!
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC MAX
!     ..
!     .. Local Scalars ..
   DOUBLE PRECISION TEMP
   INTEGER I,INFO,J,K,NROWA
   LOGICAL LSIDE,NOUNIT,UPPER
!     ..
!     .. Parameters ..
   DOUBLE PRECISION ONE,ZERO
   PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!
!     Test the input parameters.
!
   LSIDE = LSAME(SIDE,'L')
   IF (LSIDE) THEN
      NROWA = M
   ELSE
      NROWA = N
   END IF
   NOUNIT = LSAME(DIAG,'N')
   UPPER = LSAME(UPLO,'U')
!
   INFO = 0
   IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
      INFO = 1
   ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
      INFO = 2
   ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
        (.NOT.LSAME(TRANSA,'T')) .AND.&
        (.NOT.LSAME(TRANSA,'C'))) THEN
      INFO = 3
   ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
      INFO = 4
   ELSE IF (M.LT.0) THEN
      INFO = 5
   ELSE IF (N.LT.0) THEN
      INFO = 6
   ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
      INFO = 9
   ELSE IF (LDB.LT.MAX(1,M)) THEN
      INFO = 11
   END IF
   IF (INFO.NE.0) THEN
      CALL XERBLA('DTRMM ',INFO)
      RETURN
   END IF
!
!     Quick return if possible.
!
   IF (M.EQ.0 .OR. N.EQ.0) RETURN
!
!     And when  alpha.eq.zero.
!
   IF (ALPHA.EQ.ZERO) THEN
      DO 20 J = 1,N
         DO 10 I = 1,M
            B(I,J) = ZERO
10       CONTINUE
20    CONTINUE
      RETURN
   END IF
!
!     Start the operations.
!
   IF (LSIDE) THEN
      IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*A*B.
!
         IF (UPPER) THEN
            DO 50 J = 1,N
               DO 40 K = 1,M
                  IF (B(K,J).NE.ZERO) THEN
                     TEMP = ALPHA*B(K,J)
                     DO 30 I = 1,K - 1
                        B(I,J) = B(I,J) + TEMP*A(I,K)
30                   CONTINUE
                     IF (NOUNIT) TEMP = TEMP*A(K,K)
                     B(K,J) = TEMP
                  END IF
40             CONTINUE
50          CONTINUE
         ELSE
            DO 80 J = 1,N
               DO 70 K = M,1,-1
                  IF (B(K,J).NE.ZERO) THEN
                     TEMP = ALPHA*B(K,J)
                     B(K,J) = TEMP
                     IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                     DO 60 I = K + 1,M
                        B(I,J) = B(I,J) + TEMP*A(I,K)
60                   CONTINUE
                  END IF
70             CONTINUE
80          CONTINUE
         END IF
      ELSE
!
!           Form  B := alpha*A**T*B.
!
         IF (UPPER) THEN
            DO 110 J = 1,N
               DO 100 I = M,1,-1
                  TEMP = B(I,J)
                  IF (NOUNIT) TEMP = TEMP*A(I,I)
                  DO 90 K = 1,I - 1
                     TEMP = TEMP + A(K,I)*B(K,J)
90                CONTINUE
                  B(I,J) = ALPHA*TEMP
100            CONTINUE
110         CONTINUE
         ELSE
            DO 140 J = 1,N
               DO 130 I = 1,M
                  TEMP = B(I,J)
                  IF (NOUNIT) TEMP = TEMP*A(I,I)
                  DO 120 K = I + 1,M
                     TEMP = TEMP + A(K,I)*B(K,J)
120               CONTINUE
                  B(I,J) = ALPHA*TEMP
130            CONTINUE
140         CONTINUE
         END IF
      END IF
   ELSE
      IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*B*A.
!
         IF (UPPER) THEN
            DO 180 J = N,1,-1
               TEMP = ALPHA
               IF (NOUNIT) TEMP = TEMP*A(J,J)
               DO 150 I = 1,M
                  B(I,J) = TEMP*B(I,J)
150            CONTINUE
               DO 170 K = 1,J - 1
                  IF (A(K,J).NE.ZERO) THEN
                     TEMP = ALPHA*A(K,J)
                     DO 160 I = 1,M
                        B(I,J) = B(I,J) + TEMP*B(I,K)
160                  CONTINUE
                  END IF
170            CONTINUE
180         CONTINUE
         ELSE
            DO 220 J = 1,N
               TEMP = ALPHA
               IF (NOUNIT) TEMP = TEMP*A(J,J)
               DO 190 I = 1,M
                  B(I,J) = TEMP*B(I,J)
190            CONTINUE
               DO 210 K = J + 1,N
                  IF (A(K,J).NE.ZERO) THEN
                     TEMP = ALPHA*A(K,J)
                     DO 200 I = 1,M
                        B(I,J) = B(I,J) + TEMP*B(I,K)
200                  CONTINUE
                  END IF
210            CONTINUE
220         CONTINUE
         END IF
      ELSE
!
!           Form  B := alpha*B*A**T.
!
         IF (UPPER) THEN
            DO 260 K = 1,N
               DO 240 J = 1,K - 1
                  IF (A(J,K).NE.ZERO) THEN
                     TEMP = ALPHA*A(J,K)
                     DO 230 I = 1,M
                        B(I,J) = B(I,J) + TEMP*B(I,K)
230                  CONTINUE
                  END IF
240            CONTINUE
               TEMP = ALPHA
               IF (NOUNIT) TEMP = TEMP*A(K,K)
               IF (TEMP.NE.ONE) THEN
                  DO 250 I = 1,M
                     B(I,K) = TEMP*B(I,K)
250               CONTINUE
               END IF
260         CONTINUE
         ELSE
            DO 300 K = N,1,-1
               DO 280 J = K + 1,N
                  IF (A(J,K).NE.ZERO) THEN
                     TEMP = ALPHA*A(J,K)
                     DO 270 I = 1,M
                        B(I,J) = B(I,J) + TEMP*B(I,K)
270                  CONTINUE
                  END IF
280            CONTINUE
               TEMP = ALPHA
               IF (NOUNIT) TEMP = TEMP*A(K,K)
               IF (TEMP.NE.ONE) THEN
                  DO 290 I = 1,M
                     B(I,K) = TEMP*B(I,K)
290               CONTINUE
               END IF
300         CONTINUE
         END IF
      END IF
   END IF
!
   RETURN
!
!     End of DTRMM .
!
 END SUBROUTINE DTRMM
!
! SUBROUTINES included (incomplete):
! DCOPY
! DGEMM
! DGEMV
! DGETRF
! DGETRF2
! DGETRS
! DLACPY
! DLAE2
! DLAEDA
! DLAED0
! DLAED1
! DLAED2
! DLAED3
! DLAED4
! DLAED7
! DLAED8
! DLAED9
! DLAEV2
! DLAMRG
! DLARF
! DLARFG
! DLARTG
! DLASCL
! DLASET
! DLASR
! DLASRT
! DLASSQ
! DLASWP
! DLASYF
! DOPGTR
! DOPMTR
! DSCAL
! DSPEV
! DSPEVD
! DSPTRD
! DSTEDC
! DSTEQR
! DSTERF
! DSWAP
! DSYMV
! DSYR
! DSYTF2
! DSYTRF
! DSYTRI
! DTRSM
! XERBLA
!
! FUNCTIONS included
! DDOT
! DISNAN
! DLAISNAN
! DLAMCH
! DLAMC3
! DLANSP
! DLANST
! DLAPY2
! IDAMAX
! IEEECK
! ILAENV
! ILADLC
! ILADLR
! IPARMQ
! LSAME
! LSAMEN
!
END MODULE OCLABLAS


