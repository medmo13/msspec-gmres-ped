      SUBROUTINE COEF_RENORM(NDIF)
C
C  This subroutine computes the coefficients for the renormalization
C     of the multiple scattering series. These coefficients are 
C     expressed as C_REN(K) where K is the multiple scattering order.
C     REN2 is the value of the mixing (or renormalization) parameter.
C
C     NDIF is the scattering order at which the series is truncated, 
C     so that K varies from 0 to NDIF.
C
C  COMMON /RENORM/:
C
C     I_REN = 1 : renormalization in terms of G_n matrices (n : N_REN)
C           = 2 : renormalization in terms of the Sigma_n matrices
C           = 3 : renormalization in terms of the Z_n matrices
C           = 4 : renormalization in terms of the Pi_1 matrix
C
C     N_REN = n
C     
C     REN   = REN_R+IC*REN_I : omega
C
C                                        Last modified : 11 Apr 2019
C
      USE DIM_MOD
      USE C_RENORM_MOD
      USE RENORM_MOD
C
      REAL C(0:NDIF_M,0:NDIF_M)
C
      COMPLEX REN,REN2,COEF1,COEF2,ZEROC,ONEC,IC
      COMPLEX Y1(0:NDIF_M,0:NDIF_M)
C
C
      ZEROC=(0.,0.)
      ONEC=(1.,0.)
      IC=(0.,1.)
C
      REN=REN_R+IC*REN_I                                         ! omega
C
C  Initialisation of renormalization coefficients
C
      DO J=0,NDIF
        C_REN(J)=ZEROC
      ENDDO
C                                               
C  Computing the binomial coefficients C(N,K) = (N) = N! / K! (N-K)!
C                                               (K)
CCCC  2019.06.09 Aika
      c=0.0
CCCC  2019.06.09 Aika
      C(0,0)=1.
      C(1,0)=1.
      C(1,1)=1.
      DO N=2,NDIF
        C(N,0)=1.
        C(N,N)=1.
        DO K=1,N-1
          C(N,K)=C(N-1,K)+C(N-1,K-1)
        ENDDO
      ENDDO
C
      IF(I_REN.LE.3) THEN
C
C  Computing the modified renormalization parameter REN2 (g_n,s_n,zeta_n)
C
        IF(I_REN.EQ.1) THEN
C
C.....(g_n,G_n) renormalization
C
          REN2=REN**N_REN                                        ! g_n = omega^n
C
        ELSEIF(I_REN.EQ.2) THEN
C
C.....(s_{n},Sigma_n) renormalization
C
          REN2=(ONEC-REN**(N_REN+1))/(FLOAT(N_REN+1)*(ONEC-REN)) ! s_n
C
        ELSEIF(I_REN.EQ.3) THEN
C
C.....(zeta_{n},Z_n) renormalization
C
C     2019.04.29
C         REN2=(REN-ONEC)**(N_REN+1)                             ! zeta_n  
C     2019.06.09
C          REN2=-(REN-ONEC)**(N_REN+1)                             ! zeta_n  
           REN2=-(REN-ONEC)**(N_REN+1)                              ! zeta_n
C
        ENDIF
C
C     AT & MTD 2019.04.17 - summation over j ?
        DO K=0,NDIF
           c_ren(k)=zeroc
           DO J=K,NDIF
              C_REN(K)=C_REN(K)+c(j,k)*(ONEC-REN2)**(J-K)
           ENDDO
           c_ren(k)=c_ren(k)*ren2**(k+1)
        ENDDO
C        
C        DO K=0,NDIF
C          COEF1=REN2**(K+1)
C          DO J=K,NDIF
C            COEF2=(ONEC-REN2)**(J-K)
C            C_REN(J)=C_REN(J)+COEF1*COEF2*C(J,K)
C          ENDDO
C        ENDDO
C
      ELSEIF(I_REN.EQ.4) THEN
C
C     Loewdin (Pi_1) renormalization for n = 1
C
C     Notation: Y1(M,K) : [Y_1^m]_k
C     2019.06.09    
C     Notation: Y1(K,M) : [Y_1^k]_m         
C
        COEF1=ONEC-REN                                           ! (1 - omega)
        DO K=0,NDIF
          Y1(K,K)=COEF1**K
          IF(K.LT.NDIF) THEN
            DO M=K+1,NDIF
               COEF2=(REN**(M-K))*(COEF1**(2*K-M))
C 2019.04.19 AT & MTD
C            Y1(M,K)=COEF2*(C(K,M-K)+COEF1*C(K,M-K-1))
              Y1(K,M)=COEF2*(C(K,M-K)+COEF1*C(K,M-K-1))
            ENDDO
          ENDIF
        ENDDO
C
        DO K=0,NDIF
          IN=INT(K/2)
          C_REN(K)=ZEROC
          DO M=IN,K
            C_REN(K)=C_REN(K)+Y1(M,K)
          ENDDO
        ENDDO
C
      ENDIF
C
      END

C=======================================================================
C
      SUBROUTINE COEF_SIGMAPI2(NDIF)
C
C  Experimental n=2 averaged-product/Pi blend for PED path expansion.
C  M = (1-alpha) P2 + alpha (P0 + P1 + P2)/3,
C  with Pj=product(G0..Gj),
C  G0=K and Gj=(1-omega)I+omega*G(j-1).  Since M(1)=1,
C  N=(1-M)/(1-K) is a polynomial and (I-K)^-1=(I-M)^-1 N.
C  N_REN is the highest power in the finite sum: sum_{m=0}^{N_REN} M^m N.
C  Coefficients are retained through the PED path cutoff NDIF.
C
      USE DIM_MOD
      USE C_RENORM_MOD
      USE RENORM_MOD
      INTEGER D,J,IT,I,L
      REAL ALPHA
      COMPLEX REN,ONEC,ZEROC
      COMPLEX G(0:NDIF_M),GN(0:NDIF_M)
      COMPLEX P(0:NDIF_M),PN(0:NDIF_M)
      COMPLEX PS(0:NDIF_M),M(0:NDIF_M),NUM(0:NDIF_M)
      COMPLEX TERM(0:NDIF_M),TN(0:NDIF_M),SUMC(0:NDIF_M)
      COMPLEX ACC

      ONEC=(1.,0.)
      ZEROC=(0.,0.)
      REN=REN_R+(0.,1.)*REN_I
      D=NDIF
      ALPHA=0.5
      IF(I_REN.EQ.7) ALPHA=0.25
      IF(I_REN.EQ.8) ALPHA=0.75

      DO I=0,D
        G(I)=ZEROC
        GN(I)=ZEROC
        P(I)=ZEROC
        PN(I)=ZEROC
        PS(I)=ZEROC
        M(I)=ZEROC
        NUM(I)=ZEROC
        TERM(I)=ZEROC
        TN(I)=ZEROC
        SUMC(I)=ZEROC
      ENDDO
      DO I=0,NDIF
        C_REN(I)=ZEROC
      ENDDO

C  Start with G0=K and P=I; accumulate P0, P1, P2.
      IF(D.GE.1) G(1)=ONEC
      P(0)=ONEC
      DO J=0,2
        IF(J.GT.0) THEN
          DO I=0,D
            GN(I)=REN*G(I)
          ENDDO
          GN(0)=GN(0)+(ONEC-REN)
          DO I=0,D
            G(I)=GN(I)
          ENDDO
        ENDIF
        DO I=0,D
          PN(I)=ZEROC
        ENDDO
        DO I=0,D
          DO L=0,D-I
            PN(I+L)=PN(I+L)+P(I)*G(L)
          ENDDO
        ENDDO
        DO I=0,D
          P(I)=PN(I)
          PS(I)=PS(I)+P(I)
        ENDDO
      ENDDO

C  Blend product and arithmetic mean of the three partial products.
      DO I=0,D
        M(I)=(1.0-ALPHA)*P(I)+ALPHA*PS(I)/3.0
      ENDDO

C  Polynomial division (1-M)/(1-K), using cumulative coefficients.
      ACC=ONEC-M(0)
      NUM(0)=ACC
      DO I=1,D
        ACC=ACC-M(I)
        NUM(I)=ACC
      ENDDO

C  Sum M^j N through j=N_REN, truncated by the path cutoff.
      DO I=0,D
        TERM(I)=NUM(I)
        SUMC(I)=NUM(I)
      ENDDO
      DO IT=1,N_REN
        DO I=0,D
          TN(I)=ZEROC
        ENDDO
        DO I=0,D
          DO L=0,D-I
            TN(I+L)=TN(I+L)+TERM(I)*M(L)
          ENDDO
        ENDDO
        DO I=0,D
          TERM(I)=TN(I)
          SUMC(I)=SUMC(I)+TERM(I)
        ENDDO
      ENDDO
      DO I=0,D
        C_REN(I)=SUMC(I)
      ENDDO
      END
C
C=======================================================================
C
      SUBROUTINE COEF_LOEWDIN(NDIF)
C
C  This subroutine computes the coefficients for the Loewdin expansion
C     of the multiple scattering series. These coefficients are 
C     expressed as C_LOW(K) where K is the multiple scattering order.
C     REN is the value of the mixing (or renormalization) parameter.
C     NDIF is the scattering order at which the series is truncated, 
C     so that K varies from 0 to NDIF.
C
C  Corresponds to parameter I_REN = 5
C
C  Notation: X(K,N) = X_n(omega,k)
C
C
C                                        Last modified : 11 Apr 2019
C
      USE DIM_MOD
      USE C_RENORM_MOD, C_LOW => C_REN
      USE RENORM_MOD
C
      COMPLEX X(0:NDIF_M,0:NDIF_M),SUM_L,POWER
      COMPLEX REN,ZEROC,ONEC,IC
C
C
      ZEROC=(0.,0.)
      ONEC=(1.,0.)
      IC=(0.,1.)
C
      REN=REN_R+IC*REN_I                                         ! omega
C
C  Initialisation of renormalization coefficients
C
      DO J=0,NDIF
        C_LOW(J)=ZEROC
      ENDDO
C
C  Computing the X(N,K) coefficients, with K <= N
C
      POWER=ONEC/REN
      DO N=0,NDIF
        POWER=POWER*REN                                          ! omega^n
        IF(N.EQ.0) THEN 
          X(N,0)=ONEC
        ELSE
          X(N,0)=ZEROC
        ENDIF
        DO K=1,NDIF
          IF(K.GT.N) THEN
            X(N,K)=ZEROC
          ELSEIF(K.EQ.N) THEN
            X(N,K)=POWER*X(N-1,K-1)
          ELSE
            X(N,K)=X(N-1,K)*(REN-POWER) + POWER*X(N-1,K-1)
          ENDIF
        ENDDO
      ENDDO
C
C  Calculation of L_n(omega,NDIF)
C
      DO N=0,NDIF
        SUM_L=ZEROC
        DO K=N,NDIF
          SUM_L=SUM_L+X(K,N)
        ENDDO
        C_LOW(N)=SUM_L
      ENDDO

      END
