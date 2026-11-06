C
C=======================================================================
C
      SUBROUTINE COEFFICIENTS(N_APP,COEF,BETA,I_SWITCH)
C
C       Binomial coefficients     ---> ISWITCH = 1
C       Power coefficients        ---> ISWITCH = 2
C       Pochhammer coefficients   ---> ISWITCH = 3
C       Pochhammer coefficients   ---> ISWITCH = 4 (negative argument)
C
C
C  Author : D. Sebilleau
C
C                                 Last modified : 28 Sep 2013
C
C      INCLUDE 'spec.inc'
C
      USE DIM_MOD
C
      COMPLEX*16 COEF(-1:3*N_ORD_M,-1:N_ORD_M),BETA,ONEC
C
      ONEC=(1.D0,0.D0)
C
      COEF(0,0)=ONEC
      COEF(1,0)=ONEC
      COEF(1,1)=ONEC
C
      IF(I_SWITCH.EQ.1) THEN
C
        DO N=2,N_APP-1
C
          COEF(N,0)=ONEC
          COEF(N,1)=DCMPLX(N)
C
          DO K=2,N
C
            COEF(N,K)=COEF(N-1,K-1)+COEF(N-1,K)
C
          ENDDO
C
        ENDDO
C
      ELSEIF(I_SWITCH.EQ.2) THEN
C
        DO N=0,N_APP
C
          COEF(N,0)=ONEC
          COEF(N,1)=BETA+DFLOAT(N)
C
          DO K=1,N_APP
C
            COEF(N,K)=COEF(N,K-1)*(BETA+DFLOAT(N))
C
          ENDDO
C
        ENDDO
C
      ELSEIF(I_SWITCH.EQ.3) THEN
C
        DO N=0,N_APP
C
          COEF(N,0)=ONEC
C
          DO K=1,N_APP
C
            COEF(N,K)=COEF(N,K-1)*(BETA+DFLOAT(N+K-1))
C
          ENDDO
C
        ENDDO
C
      ELSEIF(I_SWITCH.EQ.4) THEN
C
        DO N=0,N_APP
C
          COEF(N,0)=ONEC
C
          DO K=1,N_APP
C
            COEF(N,K)=COEF(N,K-1)*(-BETA-DFLOAT(N+K-1))
C
          ENDDO
C
        ENDDO
C
      ENDIF
C
      RETURN
C
      END
