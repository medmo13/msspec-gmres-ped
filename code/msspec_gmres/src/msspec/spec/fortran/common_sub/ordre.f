C
C=======================================================================
C
      SUBROUTINE ORDRE(NINI,VALINI,NFIN,VALFIN)
C
C  Given a set of **real** numbers VALINI, this routine orders them and
C       suppresses the values appearing more than once. The remaining
C       values are stored in VALFIN.
C
C       VALINI(K+1).GT.VALINI(K) : decreasing order
C       VALINI(K+1).LT.VALINI(K) : increasing order
C
C
      DIMENSION VALINI(NINI),VALFIN(NINI)
C
      LOGICAL BUBBLE
C
      DATA SMALL /0.00001/
C
      DO J=1,NINI-1
         K=J
         BUBBLE=.TRUE.
150      IF(K.GE.1.AND.BUBBLE) THEN
            IF(VALINI(K+1).GT.VALINI(K)) THEN
              R1=VALINI(K)
              VALINI(K)=VALINI(K+1)
              VALINI(K+1)=R1
           ELSE
             BUBBLE=.FALSE.
           END IF
           K=K-1
           GOTO 150
         ENDIF
      ENDDO
C
      JFIN=1
      VALFIN(1)=VALINI(1)
      DO J=1,NINI-1
        IF(ABS(VALFIN(JFIN)-VALINI(J+1)).GT.SMALL) THEN
          JFIN=JFIN+1
          VALFIN(JFIN)=VALINI(J+1)
        ENDIF
      ENDDO
      NFIN=JFIN
C
      RETURN
C
      END
