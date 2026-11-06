C
C=======================================================================
C
      SUBROUTINE CHECK_CONV(J_INI,J_FIN,N_TABLE,I_CONV,N_ACC,K_ACC,
     &K_MAX,N_COEF,I_COL,ACC,RHO,SN)
C
C  This subroutine checks the convergence of the acceleration process.
C  For convergence, there must be a N_TABLE x N_TABLE square in the (n,k) table
C  of SN(n,k)that contains equal values (SMALL being the accepted error)
C
C  This way, both the column convergence (k fixed and n increasing)
C  and the diagonal convergence (n=k increasing) are considered
C
C  Important note: in some algorithms, such as the epsilon one, only
C                  odd columns are meaningful. If I_COL = 0, all
C                  columns are taken into account while for I_COL = 2
C                  only even columns are considered
C
C  Author : D. Sebilleau
C
C                                 Last modified : 27 Feb 2011
C
C      INCLUDE 'spec.inc'
C
      USE DIM_MOD
C
      COMPLEX*16 SN(-1:N_ORD_M,-1:N_ORD_M)
C
      REAL*8 ACC,REF,COMP,RHO
C
C  Step for the columns comparison and change of the size
C    of the square when only even columns should be considered
C
      IF(I_COL.EQ.0) THEN
        I_STEP=1
        N_SQUARE=N_TABLE
        N_START=N_TABLE
      ELSEIF(I_COL.EQ.2) THEN
        I_STEP=2
        N_SQUARE=N_TABLE+N_TABLE-1
        IF(MOD(N_SQUARE,2).EQ.0) THEN
          N_START=N_SQUARE
        ELSE
          N_START=N_SQUARE+1
        ENDIF
      ENDIF
C
C  Positionning the bottom right corner of the square: (N,K),
C    starting from (N_START,N_START). If only even columns
C    should be considered, N_START must be even
C
      K_FIN=J_FIN/N_COEF
      IF(K_FIN.LT.N_START) THEN
         I_CONV=0
         GOTO 10
      ENDIF
C
      DO K=N_START,K_FIN+N_START-1,I_STEP
        IF(K.GT.K_MAX) GOTO 20
        DO N=N_START,K_FIN-K,I_STEP
C
C  Checking the values (N1,K1) in the selected square
C
          REF=CDABS(SN(N,K))
          I_CONV=1
C
          DO N1=0,N_START,I_STEP
            NN=N-N1
            DO K1=0,N_START,I_STEP
C
              KK=K-K1
              COMP=CDABS(REF-SN(NN,KK))
              IF(COMP.GT.ACC) THEN
                I_CONV=0
              ENDIF
C
            ENDDO
          ENDDO
C
          IF(I_CONV.EQ.1) THEN
C
C  All values in the square are equal (within ACC):
C     Spectral radius taken as the left top corner value
C
            N_ACC=N-N_START+I_STEP
            K_ACC=K-N_START+I_STEP
            RHO=REF
            GOTO 10
          ENDIF
C
        ENDDO
  20  CONTINUE
      ENDDO
C
  10  RETURN
C
      END
