C
C=======================================================================
C
      SUBROUTINE NUMAT(NUM,NIVA,IL,IM,IN)
C
      USE OUTUNITS_MOD
      DIMENSION I(100)
C
C
      L=2*NIVA+1
      IF(L.GT.100) THEN
        WRITE(IUO1,5)
        STOP
      ENDIF
      L1=NIVA+1
C
      DO K=1,L
        IF(K.LE.L1) THEN
          I(K)=K-1
        ELSE
          I(K)=L1-K
        ENDIF
      ENDDO
C
      Q1=FLOAT(NUM)/FLOAT(L*L)
      JR1=NUM-L*L*INT(Q1+0.0001)
      JS1=INT(Q1+0.9999)
      Q2=FLOAT(JR1)/FLOAT(L)
      JS2=INT(Q2+0.9999)
      IF(JR1.EQ.0) JS2=L
      Q3=FLOAT(NUM)/FLOAT(L)
      JR3=INT(Q3+0.0001)
      JS3=NUM-L*JR3
      IF(JS3.EQ.0) JS3=L
      IL=I(JS1)
      IM=I(JS2)
      IN=I(JS3)
C
   5  FORMAT(///,'<<<<<<<<<<  INCREASE THE SIZE OF I IN','  THE NUMAT SU
     &BROUTINE  >>>>>>>>>>')
C
      RETURN
C
      END
