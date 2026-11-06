C
C=======================================================================
C
      SUBROUTINE TEST1(COUPUR,NB1,NB2,ATOME,COORD,VAL,NBAT,IRE,NBZ)
C
      USE DIM_MOD
C
      DIMENSION ATOME(3,NATCLU_M),COORD(3,NATCLU_M),VAL(NATCLU_M)
      DIMENSION IRE(NATCLU_M,2)
C
      DIST2=0.
      DO 10 I=1,3
        DIST2=DIST2+ATOME(I,NB1)*ATOME(I,NB1)
  10  CONTINUE
      DIST=SQRT(DIST2)
      V=0.0001
      IF((ATOME(3,NB1).LE.V).AND.(DIST.LE.COUPUR)) THEN
        NBAT=NBAT+1
        NB2=NB2+1
        IRE(NB1,1)=NB2
        IRE(NB1,2)=NBAT
        DO 20 I=1,3
          COORD(I,NB2)=ATOME(I,NB1)
  20    CONTINUE
        IF(NBZ.EQ.0) THEN
          NBZ=NBZ+1
          VAL(NBZ)=COORD(3,NB2)
        ELSE
          N1=0
          DO N=1,NBZ
            D=ABS(COORD(3,NB2)-VAL(N))
            IF(D.LT.0.0001) N1=N1+1
          ENDDO
          IF(N1.EQ.0) THEN
            NBZ=NBZ+1
            VAL(NBZ)=COORD(3,NB2)
          ENDIF
        ENDIF
      ENDIF
C
      RETURN
C
      END
