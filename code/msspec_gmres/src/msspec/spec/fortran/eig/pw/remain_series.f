C
C=======================================================================
C
      SUBROUTINE REMAIN_SERIES(I_GN,N_APP,SN,XN,GN)
C
C  This subroutine computes the G_N series used to describe
C     the correction on the remainder estimate
C
C      I_GN = 1 :  G_{i}(n) = G_{i-1}(n)*x_{n}
C      I_GN = 2 :  G_{i}(n) = x_{n+i-1} (G transformation)
C      I_GN = 3 :  G_{i}(n) = DELTA S_{n+i-1} (epsilon alg.)
C      I_GN = 4 :  G_{i}(n) = G_{i-1}(n)*S_{n}/x_{n}
C      I_GN = 5 :  G_{i}(n) = G_{i-1}(n)/x_{n}
C      I_GN = 6 :  G_{i}(n) = (DELTA S_{n})^i (Germain-Bonne)
C      I_GN = 7 :  G_{i}(n) = x_{n+i-2} (p process)
C
C
C  Author : D. Sebilleau
C
C                                 Last modified : 28 Feb 2013
C
C      INCLUDE 'spec.inc'
C
      USE DIM_MOD
C
      COMPLEX*16 SN(-1:N_ORD_M,-1:N_ORD_M),XN(-1:N_ORD_M)
      COMPLEX*16 GN(-1:N_ORD_M,-1:N_ORD_M,-1:N_ORD_M)
      COMPLEX*16 CNK(-1:3*N_ORD_M,-1:N_ORD_M)
      COMPLEX*16 SUM_J,MONE,ONE,ZEROC,ONEC
C
      ZEROC=(0.D0,0.D0)
      ONEC=(1.D0,0.D0)
C
C  Computing binomial coefficients
C
      IF(I_GN.EQ.6) CALL COEFFICIENTS(N_APP,CNK,ZEROC,1)
C
      DO N=0,N_APP
C
        GN(N,0,-1)=ZEROC
        GN(N,0,0)=ONEC
        MONE=-ONEC
C
        DO I=1,N_APP
C
          IF(I_GN.EQ.1) THEN
            GN(N,0,I)=GN(N,0,I-1)*XN(N)
          ELSEIF(I_GN.EQ.2) THEN
            GN(N,0,I)=XN(N+I-1)
          ELSEIF(I_GN.EQ.3) THEN
            GN(N,0,I)=SN(N+I,0)-SN(N+I-1,0)
          ELSEIF(I_GN.EQ.4) THEN
            GN(N,0,I)=GN(N,0,I-1)*SN(N,0)/XN(N)
          ELSEIF(I_GN.EQ.5) THEN
            GN(N,0,I)=GN(N,0,I-1)/XN(N)
          ELSEIF(I_GN.EQ.6) THEN
            MONE=-MONE
            SUM_J=ZEROC
            ONE=-ONEC
            DO J=0,I
              ONE=-ONE
              SUM_J=SUM_J+ONE*CNK(I,J)*SN(N+J,0)
            ENDDO
            GN(N,0,I)=MONE*SUM_J
          ELSEIF(I_GN.EQ.7) THEN
            IF(I.GT.1) THEN
              GN(N,0,I)=SN(N+I-1,0)-SN(N+I-2,0)
            ELSE
              GN(N,0,I)=SN(N,0)
            ENDIF
          ENDIF
C
        ENDDO
C
      ENDDO
C
      RETURN
C
      END
