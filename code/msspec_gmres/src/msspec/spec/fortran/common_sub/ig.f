C
C=======================================================================
C
      INTEGER FUNCTION IG(J)
C
C  This function is returns the value 1 if J is an integer
C   and 2 if it is a half-integer
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      REAL*8 J,JJ
C
      DATA SMALL /0.0001D0/
C
      JJ=ABS(J+J)
C
      LL=INT(JJ+SMALL)
C
      IF(MOD(LL,2).EQ.0) THEN
        IG=1
      ELSE
        IG=2
      ENDIF
C
      END
