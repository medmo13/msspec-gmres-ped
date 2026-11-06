C
C=======================================================================
C
      SUBROUTINE SUP_ZEROS(TL,LMAX,NE,NAT,IUO1,ITRTL)
C
C   This routine suppresses possible zeros in the TL arrays so that
C     the code runs faster because of lower values of LMAX. Actually,
C     the TL array is not modified, it is just the LMAX array that is
C     altered. This is particularly useful for energy variations or
C     for matrix inversion
C
      USE DIM_MOD
C
      COMPLEX TL_,TL(0:NT_M,4,NATM,NE_M)
C
      INTEGER LMAX(NATM,NE_M)
C
      IF(ITRTL.EQ.1) THEN
        SMALL=0.1
      ELSEIF(ITRTL.EQ.2) THEN
        SMALL=0.01
      ELSEIF(ITRTL.EQ.3) THEN
        SMALL=0.001
      ELSEIF(ITRTL.EQ.4) THEN
        SMALL=0.0001
      ELSEIF(ITRTL.EQ.5) THEN
        SMALL=0.00001
      ELSEIF(ITRTL.EQ.6) THEN
        SMALL=0.000001
      ELSEIF(ITRTL.EQ.7) THEN
        SMALL=0.0000001
      ELSEIF(ITRTL.EQ.8) THEN
        SMALL=0.00000001
      ELSE
        ITRTL=9
        SMALL=0.000000001
      ENDIF
C
      WRITE(IUO1,10)
      WRITE(IUO1,15) ITRTL
C
      DO JE=1,NE
        WRITE(IUO1,20) JE
        DO JAT=1,NAT
          NONZERO=0
          LM=LMAX(JAT,JE)
          DO L=0,LM
            TL_=TL(L,1,JAT,JE)
            IF((ABS(REAL(TL_)).GE.SMALL).OR.(ABS(AIMAG(TL_)).GE.SMALL)) 
     &THEN
               NONZERO=NONZERO+1
            ENDIF
          ENDDO
          LMAX(JAT,JE)=NONZERO-1
          WRITE(IUO1,30) JAT,LM,NONZERO-1
        ENDDO
      ENDDO
C
      WRITE(IUO1,40)
C
  10  FORMAT(//,'   ---> CHECK FOR ZEROS IN THE TL FILE TO REDUCE',' THE
     & AMOUNT OF COMPUTING :',/)
  15  FORMAT(/,'  (ONLY THE MATRIX ELEMENTS NON ZERO ','TO THE FIRST ',I
     &1,' DECIMAL DIGITS ARE KEPT)',/)
  20  FORMAT(/,15X,'ENERGY POINT No ',I3,/)
  30  FORMAT(8X,'PROTOTYPICAL ATOM No ',I5,'  INITIAL LMAX = ',I2,'   FI
     &NAL LMAX = ',I2)
  40  FORMAT(//)
C
      RETURN
C
      END
