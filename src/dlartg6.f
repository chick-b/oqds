*> \brief \b DLARTG6 generates a plane rotation with real cosine and real sine.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLARTG6 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlartg.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlartg.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlartg.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLARTG6( F, G, CS, SN, R )
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   CS, F, G, R, SN
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLARTG6 generate a plane rotation so that
*>
*>    [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*>    [ -SN  CS  ]     [ G ]     [ 0 ]
*>
*> This is a slower, more accurate version of the BLAS1 routine DROTG,
*> with the following other differences:
*>    F and G are unchanged on return.
*>    If G=0, then CS=1 and SN=0.
*>    If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*>       floating point operations (saves work in DBDSQR when
*>       there are zeros on the diagonal).
*>
*> If F exceeds G in magnitude, CS will be positive.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] F
*> \verbatim
*>          F is DOUBLE PRECISION
*>          The first component of vector to be rotated.
*> \endverbatim
*>
*> \param[in] G
*> \verbatim
*>          G is DOUBLE PRECISION
*>          The second component of vector to be rotated.
*> \endverbatim
*>
*> \param[out] CS
*> \verbatim
*>          CS is DOUBLE PRECISION
*>          The cosine of the rotation.
*> \endverbatim
*>
*> \param[out] SN
*> \verbatim
*>          SN is DOUBLE PRECISION
*>          The sine of the rotation.
*> \endverbatim
*>
*> \param[out] R
*> \verbatim
*>          R is DOUBLE PRECISION
*>          The nonzero component of the rotated vector.
*>
*>  This version has a few statements commented out for thread safety
*>  (machine parameters are computed on each entry). 10 feb 03, SJH.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLARTG6( K, L, M, N, CS, SN )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   K, L, M, N
      DOUBLE PRECISION   CS, R, SN
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
*     LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
      DOUBLE PRECISION   K1, L1, M1, N1
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
*     ..
*     .. Save statement ..
*     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
*     ..
*     .. Data statements ..
*     DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
*     IF( FIRST ) THEN
      SAFMIN = DLAMCH( 'S' )
      EPS = DLAMCH( 'E' )
      SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $     LOG( DLAMCH( 'B' ) ) / TWO )
      SAFMX2 = ONE / SAFMN2
*     FIRST = .FALSE.
*     END IF
      
      K1 = K
      M1 = M
      SCALE = MAX( K1, M1 )
      IF( SCALE.GE.SAFMX2 ) THEN
         COUNT = 0
 50      CONTINUE
         COUNT = COUNT + 1
         K1 = K1*SAFMN2
         M1 = M1*SAFMN2
         SCALE = MAX( K1, M1 )
         IF( SCALE.GE.SAFMX2 )
     $        GO TO 50
         L1 = L
         N1 = N
         DO 60 I = 1, COUNT
            L1 = L1*SAFMN2
            N1 = N1*SAFMN2
 60      CONTINUE
      ELSE IF( SCALE.LE.SAFMN2 ) THEN
         COUNT = 0
 70      CONTINUE
         COUNT = COUNT + 1
         K1 = K1*SAFMX2
         M1 = M1*SAFMX2
         SCALE = MAX( K1, M1 )
         IF( SCALE.LE.SAFMN2 )
     $        GO TO 70
         L1 = L
         N1 = N
         DO 80 I = 1, COUNT
            L1 = L1*SAFMX2
            N1 = N1*SAFMX2
 80      CONTINUE
      ELSE
         L1 = L
         N1 = N
      END IF
      
      F1=K1*L1+M1*N1
      G1=M1*L1-K1*N1
      IF( G1.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
      ELSE IF( F1.EQ.ZERO ) THEN
         CS = ZERO
         IF (G1 .LT. ZERO) THEN
            SN = -ONE
         ELSE
            SN = ONE
         ENDIF
      ELSE
         SCALE = MAX( F1, G1 )
         IF( SCALE.GE.SAFMX2 ) THEN
 10         CONTINUE
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( F1, G1 )
            IF( SCALE.GE.SAFMX2 )
     $           GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
 30         CONTINUE
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( F1, G1 )
            IF( SCALE.LE.SAFMN2 )
     $           GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
      END IF
      RETURN
*     
*     End of DLARTG6
*     
      END
