C
C=======================================================================
C
      SUBROUTINE ACC_SCAL(J_INI,J_FIN,METHOD,SN)
C
C  This subroutine computes the scalar convergence acceleration
C    for various methods
C
C  Note: let us consider that the (n,k) space pattern indicates
C        that we need S(n,k),...,S(n+p,k) in order to compute
C        S(n,k+1). We call N the maximum value of n for which
C        the S(n,0) are known. Then, to compute S(n,k+1), we
C        we need to know up to S(n+pk+p,0). This means that
C        the value of n is limited to N-pk-p.
C
C  Author : D. Sebilleau
C
C                                 Last modified : 14 Mar 2013
C
C      INCLUDE 'spec.inc'
C
      USE DIM_MOD
      USE CONVACC_MOD, L => LEVIN
C
      COMPLEX*16 SN(-1:N_ORD_M,-1:N_ORD_M),XN(-1:N_ORD_M)
      COMPLEX*16 NK(-1:3*N_ORD_M,-1:N_ORD_M)
      COMPLEX*16 SI(-1:N_ORD_M),TN(-1:N_ORD_M,-1:N_ORD_M)
      COMPLEX*16 GN(-1:N_ORD_M,-1:N_ORD_M,-1:N_ORD_M)
      COMPLEX*16 DELTA1,DELTA2,DELTA3,DELTA4,DELTA5
      COMPLEX*16 NUM,DEN,FIR
      COMPLEX*16 VAR,A,B,C,ZEROC,ONEC,XK
C
      REAL*8 EPS,EPSS,MUL
C
      CHARACTER*4 METHOD
C
      DATA EPS,EPSS  /1.D-12,1.D-150/
C
      ZEROC=(0.D0,0.D0)
      ONEC=(1.D0,0.D0)
C
      N_APP=J_FIN
      N_INI=0
C
      IF(METHOD.EQ.'AITK') THEN
C
C            The iterated Aitken and modified Aitken schemes
C
C       H. H. H. Homeier, Num. Alg. 8, 47 (1994)
C
C            I_VA = 1 : iterated Aitken            (p 49 eq 5)
C            I_VA = 2 : modified iterated Aitken   (p 49 eq 9)
C
        N_COUNT=0
C
        DO K=0,N_APP-1
C
          IF(I_VA.EQ.1) THEN
            XK=ONEC
          ELSEIF(I_VA.EQ.2) THEN
            XK=ONEC*DFLOAT(K+K+1)/DFLOAT(K+1)
          ENDIF
C
          DO N=N_INI,N_APP-2*K-2
C
            DELTA1=(SN(N+1,K)-SN(N,K))
            DELTA2=(SN(N+2,K)-SN(N+1,K))
            IF(CDABS(DELTA2).LT.EPSS) THEN
              DELTA3=ZEROC
            ELSE
              DELTA3=DELTA1*DELTA1/(DELTA2-DELTA1)
            ENDIF
            SN(N,K+1)=SN(N,K)-XK*DELTA3
            N_COUNT=N_COUNT+1
C
          ENDDO
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'RICH') THEN
C
C            The Richardson scheme
C
        N_COUNT=0
C
        DO K=0,N_APP-1
C
          DO N=N_INI,N_APP-K-2
C
            DELTA1=SN(N+1,0)-SN(N,0)
            DELTA2=SN(N+K+2,0)-SN(N+K+1,0)
            IF(CDABS(DELTA2-DELTA1).LT.EPSS) GOTO 10
            SN(N,K+1)=(SN(N+1,K)*DELTA1-SN(N,K)*DELTA2)/(DELTA1-
     &      DELTA2)
            N_COUNT=N_COUNT+1
C
  10        CONTINUE
C
          ENDDO
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'SALZ') THEN
C
C            The Salzer scheme
C
C      A. Hautot, http://physinfo.org/Acc_Conv/Acc_Conv_Part4.pdf (p 12)
C
        N_COUNT=0
C
        DO K=0,N_APP-1
C
          DO N=N_INI,N_APP-K-1
C
            DELTA1=DFLOAT(N+K+2)
            DELTA2=DFLOAT(N+1)
            IF(CDABS(DELTA2-DELTA1).LT.EPSS) GOTO 11
            SN(N,K+1)=(SN(N+1,K)*DELTA1-SN(N,K)*DELTA2)/(DELTA1-
     &      DELTA2)
            N_COUNT=N_COUNT+1
C
  11        CONTINUE
C
          ENDDO
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'GBWT') THEN
C
C            The GBW scheme
C
C     E. J. Weniger, Comp. Phys. Rep. 10, 189 (1989)
C                (eq 6.1-5 p33)
C
        N_COUNT=0
C
        CALL INTERP_POINTS(I_XN,2*N_APP,ALPHA,BETA,SN,XN)
C
        DO K=0,N_APP-1
C
          DO N=N_INI,N_APP-K-1
C
            DELTA1=XN(N)
            DELTA2=XN(N+K+1)
            IF(CDABS(DELTA1-DELTA2).LT.EPSS) GOTO 80
C
            SN(N,K+1)=(DELTA1*SN(N+1,K)-DELTA2*SN(N,K))/(DELTA1-
     &      DELTA2)
            N_COUNT=N_COUNT+1
C
  80        CONTINUE
C
          ENDDO
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'EPSI') THEN
C
C            The epsilon scheme
C
C          E. J. Weniger, Comp. Phys. Rep. 10, 189 (1989) (p 21 eq 4-2-1)
C          A. Salam, J. Comput. Appl. Math 46, 455 (1993) (p 456)
C
        N_COUNT=0
C
        CALL INTERP_POINTS(I_XN,N_APP,ALPHA,BETA,SN,XN)
C
        DO K=0,N_APP-1
C
          DO N=N_INI,N_APP-K-1
C
            DELTA1=(SN(N+1,K)-SN(N,K))
            IF(CDABS(DELTA1).LT.EPSS) GOTO 20
            IF(I_VA.EQ.1) THEN
              VAR=ONEC
            ELSEIF(I_VA.EQ.2) THEN
              VAR=XN(N+1)-XN(N)
            ELSEIF(I_VA.EQ.3) THEN
              VAR=XN(N+K+1)-XN(N+K)
            ENDIF
            SN(N,K+1)=SN(N+1,K-1)+VAR/DELTA1
            N_COUNT=N_COUNT+1
C
  20        CONTINUE
C
          ENDDO
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'EPSG') THEN
C
C            The generalized epsilon scheme
C
C             M. N. Barber and C. J. Hamer,
C    J. Austral. Math. Soc. (Series B) 23, 229 (1982)
C
C    DELTA1 = 0 implies SN(N,K+1) = 0 (calculation not performed)
C    DELTA2 = 0 implies SN(N+1,K) = SN(N,K)
C           ---> convergence achieved for (N,K)
C
        N_COUNT=0
C
        DO K=0,N_APP,2
C
          IF(I_VA.EQ.1) THEN
            SI(K)=ALPHA
            SI(K+1)=ALPHA
          ELSEIF(I_VA.EQ.2) THEN
            SI(K)=ZEROC
            SI(K+1)=-ONEC
          ENDIF
C
        ENDDO
C
        VAR=ONEC
C
        DO K=0,N_APP-1
C
          DO N=N_INI,N_APP-K-1
C
            NK(N,-1)=ZEROC
            DELTA2=(SN(N+1,K)-SN(N,K))
            IF(CDABS(DELTA2).LT.EPSS) GOTO 23
            NK(N,K)=SI(K)*NK(N,K-1)+VAR/DELTA2
            DELTA1=(NK(N,K)-NK(N-1,K))
            IF(CDABS(DELTA1).LT.EPSS) GOTO 23
            SN(N,K+1)=SN(N,K)+VAR/DELTA1
            N_COUNT=N_COUNT+1
  23        CONTINUE
C
          ENDDO
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'RHOA') THEN
C
C            The rho scheme
C
C     E. J. Weniger, Comp. Phys. Rep. 10, 189 (1989)
C
C          I_VA = 1 :  Osada's formulation (p 87 eq 11-2-1)
C          I_VA = 2 :  Standard rho        (p 34 eq 6.2-2)
C          I_VA = 3 :                      (p 37 eq 6.3-3) ???????
C
C        ALPHA is the decay parameter. It is
C        generally advised to take it as
C        an integer to ensure convergence
C
        IF(DREAL(ALPHA).LT.0.0D0) THEN
C
C  Drummond approximation to alpha :
C
C   (non valid for the Taylor expansion
C    of 1/1+x as in this case the denominator
C    is always zero)
C
          DO N=N_INI,N_APP
C
            DELTA1=SN(N+3,1)-SN(N+2,1)
            DELTA2=SN(N+2,1)-SN(N+1,1)
            DELTA3=SN(N+1,1)-SN(N,1)
C
            NUM=(DELTA2-DELTA3)*(DELTA1-DELTA2)
            DEN=(DELTA1-DELTA2)*DELTA2-(DELTA2-DELTA3)*DELTA1
            ALPHA=NUM/DEN-ONEC
C
          ENDDO
C
        ENDIF
C
        N_COUNT=0
C
        CALL INTERP_POINTS(I_XN,N_APP,ALPHA,BETA,SN,XN)
C
        DO K=0,N_APP-1
C
            DO N=N_INI,N_APP-K-1
C
              IF(I_VA.NE.3) THEN
                DEN=(SN(N+1,K)-SN(N,K))
              ELSE
                DELTA2=(SN(N+2,K)-SN(N+1,K))
                DELTA3=(SN(N+1,K)-SN(N,K))
                DEN=(XN(N+2*K)-XN(N+1))*DELTA3-(XN(N+2*K-1)-
     &          XN(N))*DELTA2
              ENDIF
C
              IF(CDABS(DEN).LT.EPS) GOTO 30
              IF((I_VA.EQ.3).AND.(N.EQ.(N_APP-K-1))) GOTO 30
C
              IF(I_VA.EQ.1) THEN
                NUM=(DFLOAT(K+1)+ALPHA)
              ELSEIF(I_VA.EQ.2) THEN
                NUM=XN(N+K+1)-XN(N)
              ELSEIF(I_VA.EQ.3) THEN
                NUM=(XN(N+2*K+1)-XN(N))*DELTA2*DELTA3
              ENDIF
C
              SN(N,K+1)=SN(N+1,K-1)+NUM/DEN
              N_COUNT=N_COUNT+1
C
  30          CONTINUE
C
            ENDDO
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'THET') THEN
C
C            The theta scheme
C
C     E. J. Weniger, Comp. Phys. Rep. 10, 189 (1989)
C     E. Calicetti et al, Phys. Rep. 446, 1 (2007)
C
C           I_VA = 1 : standard formulation     (eq 71 p 19)
C           I_VA = 2 : generalized formulation  (eq 10.1-10 p 74)
C
        KMAX2=(N_APP-2)/2
C
        N_COUNT=0
        CALL INTERP_POINTS(I_XN,N_APP,ALPHA,BETA,SN,XN)
C
        DO K=0,KMAX2
C
          IF(I_VA.EQ.3) GOTO 51
C
          DO N=N_INI,N_APP-2*K-1
C
            DEN=SN(N+1,2*K)-SN(N,2*K)
            IF(CDABS(DEN).LT.EPSS) GOTO 40
            IF(I_VA.EQ.1) THEN
              NUM=1.D0
            ELSEIF(I_VA.EQ.2) THEN
              NUM=XN(N+2*K+1)-XN(N)
            ENDIF
            SN(N,2*K+1)=SN(N+1,2*K-1)+NUM/DEN
            N_COUNT=N_COUNT+1
C
  40        CONTINUE
C
          ENDDO
C
          DO N=N_INI,N_APP-4*K-2
C
            DELTA2=SN(N+2,2*K)-SN(N+1,2*K)
            DELTA3=SN(N+1,2*K+1)-SN(N,2*K+1)
            DELTA4=SN(N+2,2*K+1)-SN(N+1,2*K+1)
C
            IF(I_VA.EQ.1) THEN
              NUM=DELTA2*DELTA4
              DEN=DELTA4-DELTA3
            ELSEIF(I_VA.EQ.2) THEN
              NUM=(XN(N+2*K+2)-XN(N))*DELTA2*DELTA4
              DEN=(XN(N+2*K+2)-XN(N+1))*DELTA3-(XN(N+2*K+1)-XN(
     &        N))*DELTA4
            ENDIF
            IF(CDABS(DEN).LT.EPSS) GOTO 50
C
            SN(N,2*K+2)=SN(N+1,2*K)+NUM/DEN
C
            N_COUNT=N_COUNT+1
C
  50        CONTINUE
C
          ENDDO
          GOTO 52
C
  51      DO N=N_INI,N_APP-4*K-3
C
            DELTA1=SN(N+1,K)-SN(N,K)
            DELTA2=SN(N+2,K)-SN(N+1,K)
            DELTA3=SN(N+3,K)-SN(N+2,K)
            DELTA4=DELTA3-DELTA2
            DELTA5=DELTA2-DELTA1
            DEN=DELTA3*DELTA5-DELTA1*DELTA4
            IF(CDABS(DEN).LT.EPSS) GOTO 53
C
            NUM=DELTA1*DELTA2*DELTA4
            SN(N,K+1)=SN(N+1,K)-NUM/DEN
C
            N_COUNT=N_COUNT+1
C
  53        CONTINUE
C
          ENDDO
C
  52      CONTINUE
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'LEGE') THEN
C
C            The Legendre-Toeplitz scheme
C
C      D. A. Smith and W. F. Ford, SIAM J. Numer. Anal. 16, 223 (1979)
C                   eq. (4.3), (4.4) p. 231
C
        N_COUNT=0
C
        DO K=0,N_APP-1
C
          DELTA1=DFLOAT(K+K+1)/DFLOAT(K+1)
          DELTA2=DFLOAT(K)/DFLOAT(K+1)
C
          DO N=N_INI,N_APP-K-1
            SN(N,K+1)=DELTA1*(2.D0*SN(N+1,K)-SN(N,K))-DELTA2*SN(
     &      N,K-1)
            N_COUNT=N_COUNT+1
C
          ENDDO
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'CHEB') THEN
C
C            The Chebyshev-Toeplitz scheme
C
C      D. A. Smith and W. F. Ford, SIAM J. Numer. Anal. 16, 223 (1979)
C                       eq. (4.4) - (4.9) p. 231
C
        SI(0)=ONEC
        SI(1)=3.D0*ONEC
C
        DO N=N_INI,N_APP-1
C
          TN(N,0)=SN(N,0)
          TN(N,1)=SN(N,0)+2.D0*TN(N+1,0)
          SN(N,0)=TN(N,0)/SI(0)
          SN(N,1)=TN(N,1)/SI(1)
C
        ENDDO
C
        N_COUNT=0
C
        DO K=1,N_APP-1
C
          SI(K+1)=6.D0*SI(K)-SI(K-1)
C
          DO N=N_INI,N_APP-K-1
C
            TN(N,K+1)=2.D0*TN(N,K)+4.D0*TN(N+1,K)-TN(N,K-1)
            SN(N,K+1)=TN(N,K+1)/SI(K+1)
            N_COUNT=N_COUNT+1
C
          ENDDO
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'OVER') THEN
C
C            The Overholt scheme
C
C          H. H. H. Homeier, J. Mol. Struc. (Theochem) 368, 81 (1996)
C                            (eq 36 p 85)
C
        N_COUNT=0
C
        DO K=0,N_APP-1
C
          DO N=N_INI,N_APP-K-2
C
            DELTA1=ONEC
            DELTA2=ONEC
            MUL=1.D0
C
            DO I=1,K+1
C
              DELTA1=DELTA1*(SN(N+K+1,0)-SN(N+K,0))
              DELTA2=DELTA2*(SN(N+K+2,0)-SN(N+K+1,0))
              MUL=MUL*EPS
C
            ENDDO
C
            IF(CDABS(DELTA1-DELTA2).LT.MUL) GOTO 60
C
            SN(N,K+1)=(DELTA1*SN(N+1,K)-DELTA2*SN(N,K))/(DELTA1-
     &      DELTA2)
            N_COUNT=N_COUNT+1
C
  60        CONTINUE
C
          ENDDO
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'DURB') THEN
C
C            The Durbin scheme
C
C        C. Brezinski and M. Redivo Zaglia, Comput. Appl. Math. 26, 171 (2007)
C                 (eq 25 p 185)
C
        N_COUNT=0
C
        DO K=0,N_APP-1
C
          DO N=N_INI,N_APP-4*K-4
C
            DELTA1=3.D0*SN(N+2,K)*SN(N+2,K)-4.D0*SN(N+1,K)*SN(N+
     &      3,K)+SN(N,K)*SN(N+4,K)
            DELTA2=6.D0*SN(N+2,K)-4.D0*(SN(N+1,K)+SN(N+3,K))+SN(
     &      N,K)+SN(N+4,K)
            IF(CDABS(DELTA2).LT.EPSS) GOTO 70
            SN(N,K+1)=DELTA1/DELTA2
            N_COUNT=N_COUNT+1
C
  70        CONTINUE
C
          ENDDO
C
        ENDDO
C
      ELSEIF(METHOD(2:4).EQ.'LEV') THEN
C
C            The Levin schemes
C
        CALL LEVIN(ALPHA,BETA,SN,I_WN,I_XN,N_APP,N_INI,L,N_COUNT,
     &  METHOD(1:1))
C
      ELSEIF(METHOD.EQ.'EULE') THEN
C
C            The generalized Euler scheme
C
C      D. A. Smith and W. F. Ford, SIAM J. Numer. Anal. 16, 223 (1979)
C
        N_COUNT=0
C
        DO K=0,N_APP-1
C
          DO N=N_INI,N_APP-K-1
C
            SN(N,K+1)=(SN(N+1,K)-BETA*SN(N,K))/(ONEC-BETA)
            N_COUNT=N_COUNT+1
C
          ENDDO
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'VARI') THEN
C
C            Various schemes
C
C          I_VA = 1 : Weniger lambda scheme (p 87 eq 11-2-1)
C          I_VA = 2 : Weniger sigma scheme  (p 87 eq 11-2-2)
C          I_VA = 3 : Weniger mu scheme     (p 87 eq 11-2-3)
C          I_VA = 4 : iterated rho scheme
C          I_VA = 5 : Bjorstad, Dahlquist and Grosse scheme (p 77 eq 6-3-4)
C
C          E. J. Weniger, Comp. Phys. Rep. 10, 189 (1989)
C
        N_COUNT=0
C
        CALL INTERP_POINTS(I_XN,N_APP,ALPHA,BETA,SN,XN)
C
        DO K=0,N_APP-1
C
          DO N=N_INI,N_APP-2*K-2
C
            IF(I_VA.EQ.1) THEN
              A=BETA+DFLOAT(N)
              B=A+ONEC
              C=A
            ELSEIF(I_VA.EQ.2) THEN
              A=ALPHA+DFLOAT(N+K)
              B=A+ONEC
              C=A
            ELSEIF(I_VA.EQ.3) THEN
              A=ALPHA+DFLOAT(N-K)
              B=A+ONEC
              C=A
            ELSEIF(I_VA.EQ.4) THEN
              A=XN(N+2*K+2)-XN(N)
              B=XN(N+2*K+1)-XN(N)
              C=XN(N+2*K+2)-XN(N+1)
            ELSEIF(I_VA.EQ.5) THEN
              A=(DFLOAT(2*K+1)+ALPHA)/(DFLOAT(2*K)+ALPHA)
              B=ONEC
              C=ONEC
            ENDIF
C
            DELTA1=SN(N+1,K)-SN(N,K)
            DELTA2=SN(N+2,K)-SN(N+1,K)
            DEN=B*DELTA2-C*DELTA1
            IF(CDABS(DEN).LT.EPSS) GOTO 81
            SN(N,K+1)=SN(N+1,K)-A*DELTA1*DELTA2/DEN
            N_COUNT=N_COUNT+1
C
  81        CONTINUE
C
          ENDDO
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'ITHE') THEN
C
C               Iterated theta schemes
C
C       R. Thukral, Appl. Math. Comput. 187, 1502 (2007) <--- Lubkin
C       E. J. Weniger, Comp. Phys. Rep. 10, 189 (1989)
C
C         I_VA = 1 : Lubkin transform = iterated theta transform (p 79 eq 10-3-6)
C         I_VA = 2 : iterated theta transform                    (p 84 eq 11-1-5)
C         I_VA = 3 : 2nd modified iterated Aitken transform      (p 85 eq 11-1-12)
C
        N_COUNT=0
C
        DO K=0,N_APP-1
C
          DO N=N_INI,N_APP-3*K-3
C
            DELTA1=SN(N+1,K)-SN(N,K)
            DELTA2=SN(N+2,K)-SN(N+1,K)
            DELTA3=SN(N+3,K)-SN(N+2,K)
            IF(I_VA.EQ.1) THEN
              NUM=-DELTA1*DELTA2*(DELTA3-DELTA2)
              DEN=DELTA3*(DELTA2-DELTA1)-DELTA1*(DELTA3-DELTA2)
              FIR=SN(N+1,K)
            ELSEIF(I_VA.EQ.2) THEN
              NUM=DELTA1*DELTA1*DELTA1*(DELTA3-DELTA2)
              DEN=DELTA1*DELTA1*(DELTA3-DELTA2)-DELTA2*DELTA2*(
     &        DELTA2-DELTA1)
              FIR=SN(N+1,K)
            ELSEIF(I_VA.EQ.3) THEN
              NUM=DELTA2*DELTA2*DELTA3*(DELTA3-DELTA2)
              DEN=DELTA2*DELTA2*(DELTA3-DELTA2)-DELTA3*DELTA3*(
     &        DELTA2-DELTA1)
              FIR=SN(N+2,K)
            ENDIF
C
            IF(CDABS(DEN).LT.EPSS) GOTO 90
            SN(N,K+1)=FIR+NUM/DEN
            N_COUNT=N_COUNT+1
C
  90        CONTINUE
C
          ENDDO
C
        ENDDO
C
      ELSEIF(METHOD.EQ.'EALG') THEN
C
C            E algorithm
C
        N_COUNT=0
C
        CALL INTERP_POINTS(I_XN,2*N_APP,ALPHA,BETA,SN,XN)
        CALL REMAIN_SERIES(I_GN,N_APP,SN,XN,GN)
C
C  Computing the GN(N,K,I)
C
        DO K=1,N_APP
C
          DO N=N_INI,N_APP
C
            DEN=GN(N+1,K-1,K)-GN(N,K-1,K)
            IF(CDABS(DEN).LT.EPSS) GOTO 91
            DELTA1=GN(N,K-1,K)/DEN
C
            DO I=0,N_APP
C
              NUM=GN(N+1,K-1,I)-GN(N,K-1,I)
              GN(N,K,I)=GN(N,K-1,I)-NUM*DELTA1
C
            ENDDO
C
  91        CONTINUE
C
          ENDDO
C
        ENDDO
C
C  Computing SN(N,K)
C
        DO K=1,N_APP
C
          DO N=N_INI,N_APP-K
C
            DEN=GN(N+1,K-1,K)-GN(N,K-1,K)
            IF(CDABS(DEN).LT.EPSS) GOTO 92
            DELTA1=GN(N,K-1,K)/DEN
            SN(N,K)=SN(N,K-1)-(SN(N+1,K-1)-SN(N,K-1))*DELTA1
C
  92        CONTINUE
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
