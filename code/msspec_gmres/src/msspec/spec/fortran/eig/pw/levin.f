C
C=======================================================================
C
      SUBROUTINE LEVIN(ALPHA,BETA,SN,I_WN,I_XN,N_APP,N_INI,L,N_COUNT,
     &OMEGA)
C
C  This subroutine computes Levin-type transforms of
C    a given sequence SN(N,1)
C
C     E. J. Weniger, Comp. Phys. Rep. 10, 189 (1989)
C
C  Author : D. Sebilleau
C
C                                 Last modified : 28 Feb 2013
C
C      INCLUDE 'spec.inc'
C
      USE DIM_MOD
C
      COMPLEX*16 SN(-1:N_ORD_M,-1:N_ORD_M),AN(-1:N_ORD_M),WN(-1:
     &N_ORD_M)
      COMPLEX*16 COEF(-1:3*N_ORD_M,-1:N_ORD_M),ZEROC,ONEC
      COMPLEX*16 COEL(-1:3*N_ORD_M,-1:N_ORD_M),XN(-1:N_ORD_M)
      COMPLEX*16 DEN(-1:N_ORD_M,-1:N_ORD_M),NUM(-1:N_ORD_M,-1:N_ORD_M)
      COMPLEX*16 X,XK,XXK,ALPHA,BETA,INC,SUM_1,SUM_2
C
      CHARACTER*1 OMEGA
C
      ZEROC=(0.D0,0.D0)
      ONEC=(1.D0,0.D0)
C
C  Recovery of the original series elements
C
      DO N=N_INI,N_APP
C
        AN(N)=SN(N,0)-SN(N-1,0)
C
      ENDDO
C
C
C  Choice of the remainder estimate omega_n
C
C  Reference : H. H. H. Homeier, J. Comput. App. Math. 122, 81 (2000)
C
C
      IF(OMEGA.EQ.'U') THEN
C
        MAX_N=N_APP
C
        DO N=N_INI,MAX_N
          WN(N)=AN(N)*(DFLOAT(N)+BETA)
        ENDDO
C
      ELSEIF(OMEGA.EQ.'V') THEN
C
        MAX_N=N_APP-1
C
        DO N=N_INI,MAX_N
          WN(N)=AN(N)*AN(N+1)/(AN(N)-AN(N+1))
        ENDDO
C
      ELSEIF(OMEGA.EQ.'T') THEN
C
        MAX_N=N_APP
C
        DO N=N_INI,MAX_N
          WN(N)=AN(N)
        ENDDO
C
      ELSEIF(OMEGA.EQ.'D') THEN
C
        MAX_N=N_APP-1
C
        DO N=N_INI,MAX_N
          WN(N)=AN(N+1)
        ENDDO
C
      ELSEIF(OMEGA.EQ.'E') THEN
C
        MAX_N=N_APP-2
C
        DO N=N_INI,MAX_N
          WN(N)=AN(N)*AN(N+2)/AN(N+1)
        ENDDO
C
      ENDIF
C
C  Choice of the Levin-type transform
C
C     I_WN = 1 ---> L-type                     (eq 7.2-8 p 41) (eq 3-13b p 10)
C     I_WN = 2 ---> S-type                     (eq 8.3-7 p 58) (eq 3-14b p 10)
C     I_WN = 3 ---> M-type  !! BETA > K-1      (eq 9.3-7 p 66) (eq 3-15b p 10)
C     I_WN = 4 ---> C-type                                     (eq 3-16c p 11)
C     I_WN = 5 ---> Drummond-type              (eq 9.5-5 p 70)
C     I_WN = 6 ---> R-type                     (eq 7.14-12 p 47)
C
C  References : E. J. Weniger, J. Math. Phys. 45 1209 (2004)
C               E. J. Weniger, Comp. Phys. Rep. 10 189 (1989)
C
C
C  Check for the possibility of L extension (the standard case is L = 0)
C
      IF(L.GT.0) THEN
        IF(I_WN.EQ.1) THEN
          CALL COEFFICIENTS(MAX_N,COEL,BETA,2)
          DO N=N_INI,MAX_N
            WN(N)=WN(N)*COEL(N,L)
          ENDDO
        ELSEIF(I_WN.EQ.2) THEN
          CALL COEFFICIENTS(N_APP,COEL,BETA,3)
          DO N=N_INI,MAX_N
            WN(N)=WN(N)*COEL(N,L)
          ENDDO
        ELSEIF(I_WN.EQ.3) THEN
          CALL COEFFICIENTS(N_APP,COEL,BETA,4)
          DO N=N_INI,MAX_N
            WN(N)=WN(N)*COEL(N,L)
          ENDDO
        ENDIF
      ENDIF
C
      DO N=N_INI,MAX_N
C
        X=BETA+DFLOAT(N)
        COEF(N,0)=ONEC
C
        DO K=1,N_APP
C
          XK=DFLOAT(K)
          XXK=X+XK
C
          IF(I_WN.EQ.1) THEN
            INC=XXK/(XXK+ONEC)
            COEF(N,K)=X*(INC**K)/XXK
          ELSEIF(I_WN.EQ.2) THEN
            COEF(N,K)=(XXK-ONEC)*XXK/((XXK+XK-ONEC)*(XXK+XK))
          ELSEIF(I_WN.EQ.3) THEN
            COEF(N,K)=(X-XK+ONEC)/(XXK+ONEC)
          ELSEIF(I_WN.EQ.4) THEN
            SUM_1=ZEROC
            SUM_2=ZEROC
            DO J=0,K
              SUM_1=SUM_1+ALPHA*XXK
              SUM_2=SUM_2+ALPHA*(XXK+ONEC)
            ENDDO
            COEF(N,K)=(ALPHA*X+XK-ONEC)*SUM_1/(ALPHA*XXK*SUM_2)
          ELSEIF(I_WN.EQ.5) THEN
            COEF(N,K)=ONEC
          ELSEIF(I_WN.EQ.6) THEN
            CALL INTERP_POINTS(I_XN,N_APP,ALPHA,BETA,SN,XN)
          ENDIF
C
        ENDDO
C
      ENDDO
C
C  Computation of the recurrence relation for numerator
C                and denominator
C
      DO N=N_INI,MAX_N
C
C  Starting values
C
        NUM(N,0)=SN(N,0)/WN(N)
        DEN(N,0)=ONEC/WN(N)
        SN(N,0)=NUM(N,0)/DEN(N,0)
C
      ENDDO
C
      N_COUNT=0
C
      DO K=0,N_APP-1
C
        DO N=N_INI,MAX_N-K-1
C
          IF(I_WN.NE.6) THEN
            NUM(N,K+1)=NUM(N+1,K)-COEF(N,K)*NUM(N,K)
            DEN(N,K+1)=DEN(N+1,K)-COEF(N,K)*DEN(N,K)
          ELSE
            NUM(N,K+1)=(NUM(N+1,K)-NUM(N,K))/(XN(N+K+1)-XN(N))
            DEN(N,K+1)=(DEN(N+1,K)-DEN(N,K))/(XN(N+K+1)-XN(N))
          ENDIF
          SN(N,K+1)=NUM(N,K+1)/DEN(N,K+1)
          N_COUNT=N_COUNT+1
C
        ENDDO
C
      ENDDO
C
      RETURN
C
      END
