C
C=======================================================================
C
      SUBROUTINE ORDRE2(NINI,VALINI,NFIN,VALFIN)
C
C  Given a set of **integer** numbers VALINI, this routine orders them
C       and suppresses the values appearing more than once. The remaining
C       values are stored in VALFIN.
C
C       VALINI(K+1).GT.VALINI(K) : decreasing order
C       VALINI(K+1).LT.VALINI(K) : increasing order
C
C
C
      INTEGER VALINI(NINI),VALFIN(NINI),R1
C
      LOGICAL BUBBLE
C
      DO J=1,NINI-1
         K=J
         BUBBLE=.TRUE.
150      IF(K.GE.1.AND.BUBBLE) THEN
            IF(VALINI(K+1).LT.VALINI(K)) THEN
              R1=VALINI(K)
              VALINI(K)=VALINI(K+1)
              VALINI(K+1)=R1
           ELSE
             BUBBLE=.FALSE.
           ENDIF
           K=K-1
           GOTO 150
         ENDIF
      ENDDO
C
      JFIN=1
      VALFIN(1)=VALINI(1)
      DO J=1,NINI-1
        IF(ABS(VALFIN(JFIN)-VALINI(J+1)).GT.0) THEN
          JFIN=JFIN+1
          VALFIN(JFIN)=VALINI(J+1)
        ENDIF
      ENDDO
      NFIN=JFIN
C
      RETURN
C
      END
