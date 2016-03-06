      SUBROUTINE DOQDS3(N0,A,B,WORK,INFO)

      IMPLICIT NONE

      INTEGER N0, INFO, ISUB, MAXITER, SIT
      DOUBLE PRECISION A(*), B(*), WORK(*)
      INTEGER N, M, I, J, K
      INTEGER INDRV1, INDRV2
      DOUBLE PRECISION TMP1, TMP2, TMP3, TMP4
      DOUBLE PRECISION TAU, TAU1, TAU2, TAU3
      DOUBLE PRECISION SIGMA, SIGMA2, DESIG, T, DESIG0
      DOUBLE PRECISION C1, S1, SMIN, EPS, TOL

      DOUBLE PRECISION ONE, ZERO, HALF, TWO
      INTEGER MAXREC
      PARAMETER (ONE = 1.0D0, ZERO = 0.0D0, HALF = 0.5D0, TWO = 2.0D0)
      DOUBLE PRECISION HUNDRD
      PARAMETER (HUNDRD = 100.0D0)
      PARAMETER (MAXREC = 2000000000)
*
      EXTERNAL DLAMCH
      DOUBLE PRECISION DLAMCH
      EXTERNAL DFMA0
      DOUBLE PRECISION DFMA0
*
      A(1:N0)=ABS(A(1:N0))
      B(1:N0-1)=ABS(B(1:N0-1))
*
      INDRV1 = 0
      INDRV2 = INDRV1+N0
*
      EPS = DLAMCH( 'Precision' )
      TOL = HUNDRD*EPS
*
      M = 1
      N = N0
*
      IF (N-M+1 .GE. 3) THEN

         IF (A(M) .LT. A(N)) THEN
            K=(N-M+1)/2+M-1
            DO J=M,K
               TMP1=A(J)
               A(J)=A(N+M-J)
               A(N+M-J)=TMP1
               TMP1=B(J)
               B(J)=B(N-1+M-J)
               B(N-1+M-J)=TMP1
            ENDDO
         END IF
*
         TMP1 = A(M)
         DO J = M, N-3
            IF (TMP1 .EQ. B(J)+TMP1) THEN
               B(J) = -ZERO
               WORK(INDRV2+J) = -ZERO
               M = J+1
               TMP1 = A(J+1)
            ELSE
               TMP1 = A(J+1)*(TMP1/(TMP1+B(J)))
            ENDIF
         ENDDO

         IF (TMP1 .EQ. B(N-2)+TMP1) THEN
            B(N-2) = ZERO
            TMP1 = A(N-1)
         ELSE
            TMP1 = A(N-1)*(TMP1/(TMP1+B(N-2)))
         ENDIF

         IF (TMP1 .EQ. B(N-1)+TMP1) THEN
            B(N-1) = ZERO
         ENDIF
         
      ENDIF
*
      B(N) = ZERO
      WORK(INDRV2+N) = ZERO
*
 3000 SIGMA = -B(N)
      SIGMA2 = TOL*SIGMA
      DESIG = -WORK(INDRV2+N)

      MAXITER = 30*(N-M+1)
      DO I = 1, MAXITER
*
 15      IF (N-M+1 .EQ. 1) THEN
            
            CALL DLARTG7(SIGMA,DESIG,A(N),A(N),DESIG0)
            
            GO TO 700
            
         ENDIF
         
         IF (N-M+1 .EQ. 2) THEN
            
            CALL DLAS2(A(N-1), B(N-1), A(N), A(N), A(N-1))
            
            CALL DLARTG7(SIGMA,DESIG,A(N-1),A(N-1),DESIG0)
            CALL DLARTG7(SIGMA,DESIG,A(N),A(N),DESIG0)
            
            GO TO 700
         ENDIF

         IF (B(N-1) .LE. SIGMA2) THEN
            
            CALL DLARTG7(SIGMA,DESIG,A(N),A(N),DESIG0)
            
            N = N - 1
            GO TO 15
            
         ENDIF
         
         IF (B(N-2) .LE. SIGMA2) THEN
            
            CALL DLAS2(A(N-1), B(N-1), A(N), A(N), A(N-1))
            
            CALL DLARTG7(SIGMA,DESIG,A(N-1),A(N-1),DESIG0)
            CALL DLARTG7(SIGMA,DESIG,A(N),A(N),DESIG0)
            
            N = N - 2
            GO TO 15

         ENDIF

         IF (A(M) .LT. A(N)) THEN
            K=(N-M+1)/2+M-1
            DO J=M,K
               TMP1=A(J)
               A(J)=A(N+M-J)
               A(N+M-J)=TMP1
               TMP1=B(J)
               B(J)=B(N-1+M-J)
               B(N-1+M-J)=TMP1
            ENDDO
         END IF
*     
         CALL DLAS2(A(N-1), B(N-1), A(N), TMP2, TMP3)

         TAU2 = MIN(TMP2,A(N))
         IF (TAU2 .EQ. ZERO) GO TO 350

         CALL DLARTG7(SIGMA,DESIG,TAU2,T,DESIG0)

         IF (TMP3 .GT. SQRT(TWO)*TMP2) THEN
            SIT = 1
         ELSE
            SIT = 0
         ENDIF
         IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 350
*     
         TAU1 = MINVAL(A(M:N-1))
         IF (TAU2 .GE. TAU1) THEN
            IF (TAU1 .EQ. ZERO) GO TO 350
            CALL DLARTG7(SIGMA,DESIG,TAU1,T,DESIG0)
            IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 350
            GO TO 160
         ENDIF
*     
         TAU = TAU2
         TMP4 = A(M)-TAU
         IF (TMP4 .LE. ZERO) THEN
            GO TO 160
         ELSE
            TMP3 = SQRT(TMP4)*SQRT(A(M)+TAU)
         ENDIF
         DO J = M, N-2
            CALL DLARTG4(TMP3,B(J),C1)
            TMP4 = DFMA0(C1,A(J+1),-TAU)
            IF (TMP4 .LE. ZERO) THEN
               GO TO 160
            ELSE
               TMP3 = SQRT(TMP4)*SQRT(DFMA0(C1,A(J+1),TAU))
            ENDIF
         ENDDO
         CALL DLARTG4(TMP3,B(N-1),C1)
         TMP4 = DFMA0(C1,A(N),-TAU)
         IF (TMP4 .LT. ZERO) THEN
            DO K=0,MAXREC
               TAU3 = MAX(ONE-DBLE(K)*EPS,ZERO)*(C1*A(N))
               IF (TAU3 .LT. TAU) EXIT
            ENDDO
            TAU = TAU3
         ENDIF
         IF (TAU .GT. ZERO) GO TO 125
*     
 160     WORK(INDRV1+N) = ONE/A(N)
         TMP3 = WORK(INDRV1+N)
         DO J = N-1,M,-1
            WORK(INDRV1+J) = (ONE+
     $           B(J)*WORK(INDRV1+J+1))/A(J)
            TMP3 = MAX(TMP3,WORK(INDRV1+J))
         ENDDO
         WORK(INDRV1+M) = (WORK(INDRV1+M)/TMP3)/A(M)
         TMP1 = WORK(INDRV1+M)
         DO J = M+1,N
            WORK(INDRV1+J) = (WORK(INDRV1+J)/TMP3+
     $           B(J-1)*WORK(INDRV1+J-1))/A(J)
            TMP1 = MAX(TMP1,WORK(INDRV1+J))
         ENDDO
         TAU = MAX(ZERO,(ONE/SQRT(TMP1))/SQRT(TMP3))
         
         WORK(INDRV1+N) = WORK(INDRV1+N)/TMP1
         WORK(INDRV2+N) = WORK(INDRV1+N)/A(N)
         TMP3 = WORK(INDRV2+N)
         DO J = N-1,M,-1
            WORK(INDRV1+J) = WORK(INDRV1+J)/TMP1
            WORK(INDRV2+J) = (WORK(INDRV1+J)+
     $           B(J)*WORK(INDRV2+J+1))/A(J)
            TMP3 = MAX(TMP3,WORK(INDRV2+J))
         ENDDO
         TMP2 = (WORK(INDRV2+M)/TMP3)/A(M)
         TMP1 = WORK(INDRV1+M)/TMP2
         DO J = M+1,N
            TMP2 = (WORK(INDRV2+J)/TMP3+
     $           B(J-1)*TMP2)/A(J)
            TMP1 = MIN(TMP1,WORK(INDRV1+J)/TMP2)
         ENDDO
         TAU = MAX(TAU,SQRT(TMP1)/SQRT(TMP3))
         
         IF (TAU .GT. ZERO) GO TO 125
         
         TAU=A(N)-HALF*B(N-1)
         IF (TAU .LE. ZERO) GO TO 350
         DO J=N-1,M+1,-1
            TMP2=A(J)-HALF*(B(J)+B(J-1))
            IF (TMP2 .LE. ZERO) GO TO 350
            TAU=MIN(TAU,TMP2)
         ENDDO
         TMP2=A(M)-HALF*B(M)
         IF (TMP2 .LE. ZERO) GO TO 350
         TAU=MIN(TAU,TMP2)
*     
 125     IF (TAU .EQ. ZERO) GO TO 350
         CALL DLARTG7(SIGMA,DESIG,TAU,T,DESIG0)
         IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 350
*
         TMP4 = A(M)-TAU
         IF (TMP4 .LE. ZERO) THEN
            DO K=0,MAXREC
               TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*A(M)
               IF (TAU2 .LT. TAU) EXIT
            ENDDO
            TAU = MAX(HALF*TAU,TAU2)
            GO TO 125
         ELSE
            TMP3 = SQRT(TMP4)*SQRT(A(M)+TAU)
         ENDIF
         DO J = M, N-2
            CALL DLARTG2(TMP3,B(J),C1,S1,WORK(INDRV1+J))
            WORK(INDRV2+J) = S1*A(J+1)
            TMP4 = DFMA0(C1,A(J+1),-TAU)
            IF (TMP4 .LE. ZERO) THEN
               DO K=0,MAXREC
                  TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*(C1*A(J+1))
                  IF (TAU2 .LT. TAU) EXIT
               ENDDO
               TAU = MAX(HALF*TAU,TAU2)
               GO TO 125
            ELSE
               TMP3 = SQRT(TMP4)*SQRT(DFMA0(C1,A(J+1),TAU))
            ENDIF
         ENDDO
         CALL DLARTG2(TMP3,B(N-1),C1,S1,WORK(INDRV1+N-1))
         WORK(INDRV2+N-1) = S1*A(N)
         TMP4 = DFMA0(C1,A(N),-TAU)
         IF (TMP4 .LT. ZERO) THEN
            DO K=0,MAXREC
               TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*(C1*A(N))
               IF (TAU2 .LT. TAU) EXIT
            ENDDO
            TAU = MAX(HALF*TAU,TAU2)
            GO TO 125
         ELSE
            TMP3 = SQRT(TMP4)*SQRT(DFMA0(C1,A(N),TAU))
         ENDIF
         WORK(INDRV1+N) = TMP3
*
         SIGMA = T
         SIGMA2 = TOL*SIGMA
         DESIG = DESIG0
*
         CALL DCOPY(N-M+1,WORK(INDRV1+M),1,A(M),1)
         CALL DCOPY(N-M,WORK(INDRV2+M),1,B(M),1)
*
         GO TO 400
*     
 350     TMP1 = A(M)
         DO J = M, N-1
            CALL DLARTG2(TMP1,B(J),C1,S1,A(J))
            B(J) = S1*A(J+1)
            TMP1 = C1*A(J+1)
         ENDDO
         A(N) = TMP1
*
 400     TMP1 = A(M)
         DO J = M, N-3
            IF (B(J) .LE. SIGMA2 .OR. TMP1 .EQ. B(J)+TMP1) THEN
               B(J) = -SIGMA
               WORK(INDRV2+J) = -DESIG
               M = J+1
               TMP1 = A(J+1)
            ELSE
               TMP1 = A(J+1)*(TMP1/(TMP1+B(J)))
            ENDIF
         ENDDO

         IF (B(N-2) .LE. SIGMA2 .OR. TMP1 .EQ. B(N-2)+TMP1) THEN
            B(N-2) = ZERO
            TMP1 = A(N-1)
         ELSE
            TMP1 = A(N-1)*(TMP1/(TMP1+B(N-2)))
         ENDIF

         IF (B(N-1) .LE. SIGMA2 .OR. TMP1 .EQ. B(N-1)+TMP1) THEN
            B(N-1) = ZERO
         ENDIF
*     
      ENDDO
      INFO = 2
      RETURN
         
 700  N = M-1
      IF (N .LE. 0) THEN
         GO TO 4000
      ENDIF
*
      DO J = M-2,1,-1
         IF (B(J) .LE. ZERO) THEN
            M = J+1
            GO TO 3000
         ENDIF
      ENDDO
      M = 1
      GO TO 3000
*
 4000 DO I = 1, N0 - 1
*
*     Scan for smallest D(I)
*     
         ISUB = 1
         SMIN = A( 1 )
         DO J = 2, N0 + 1 - I
            IF( A( J ).LE.SMIN ) THEN
               ISUB = J
               SMIN = A( J )
            END IF
         ENDDO
         IF( ISUB.NE.N0+1-I ) THEN
*     
*     Swap singular values and vectors
*     
            A( ISUB ) = A( N0+1-I )
            A( N0+1-I ) = SMIN
         END IF
      ENDDO
*
      INFO = 0
      RETURN
*
      END
