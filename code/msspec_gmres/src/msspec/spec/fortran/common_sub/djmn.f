C
C=======================================================================
C
      SUBROUTINE DJMN(RBETA,R,LMAX)
C
C  This routine calculates Wigner rotation matrices R^{L}_{M1 M2} up to
C             order LMAX, following Messiah's convention.
C                   They are stored as R(M2,M1,L).
C
C                                         Last modified : 20 Oct 2006
C
      USE DIM_MOD
C
      USE COEFRLM_MOD
      USE EXPROT_MOD
C
      INTEGER EPS0
C
      DIMENSION R(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1)
C
      DATA SMALL,SQR2 /0.001,1.4142136/
C
      C=COS(RBETA)*0.5
      S=SIN(RBETA)*0.5
      CC=C+C
      CMUL=-1.
      IF(ABS(S).LT.SMALL) THEN
        IF(C.GT.0.) EPS0=1
        IF(C.LT.0.) EPS0=-1
        DO L=0,LMAX
          DO M1=-L,L
            DO M2=-L,L
              IF(M1.NE.M2*EPS0) THEN
                R(M2,M1,L)=0.
              ELSE
                IF(EPS0.EQ.1) THEN
                  R(M2,M1,L)=1.
                ELSE
                  IF(MOD(L+M1,2).EQ.0) THEN
                    R(M2,M1,L)=1.
                  ELSE
                    R(M2,M1,L)=-1.
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        S1=S*SQR2
        C1=0.5+C
        R(0,0,0)=1.0
        R(-1,-1,1)=C1
        R(0,-1,1)=S1
        R(1,-1,1)=1.-C1
        R(-1,0,1)=-S1
        R(0,0,1)=CC
        R(1,0,1)=S1
        R(-1,1,1)=1.-C1
        R(0,1,1)=-S1
        R(1,1,1)=C1
C
        PRODL=-S
        COEF=-S/C1
        CL=-1.
        DO L=2,LMAX
          CL=-CL
          L1=L-1
          FLL1=CC*FLOAT(L+L1)
          FLL2=1./(FLOAT(L*L1)*CC)
          PRODL=-PRODL*S
C
C  Case M = 0
C
          R_1=EXPR(0,L)*PRODL
          R(-L,0,L)=R_1
C
          R(L,0,L)=R_1*CL
          R(0,-L,L)=R_1*CL
C
          R(0,L,L)=R_1
C
          CM2=CL
          DO M2=-L1,-1
            CM2=CM2*CMUL
            CF1=CF(L1,0,-M2)/FLL1
            CF2=FLL1/CF(L,0,-M2)
            IF(-M2.LT.L1) THEN
              R_A=CF2*(R(M2,0,L1)-R(M2,0,L-2)*CF1)
            ELSE
              R_A=CF2*R(M2,0,L1)
            ENDIF
C
            R(M2,0,L)=R_A
C
            R(-M2,0,L)=R_A*CM2
            R(0,M2,L)=R_A*CM2
C
            R(0,-M2,L)=R_A
C
          ENDDO
C
          R(0,0,L)=FLL1*R(0,0,L1)/CF(L,0,0)-R(0,0,L-2)*CF(L1,0,0)/CF(L,0
     &,0)
C
C  Case M > 0
C
          PRODM=1.
          CM=CL
          FLLM=0.
          DO M=1,L1
            CM=-CM
            PRODM=PRODM*COEF
            FLLM=FLLM+FLL2
C
            R_1=EXPR(M,L)*PRODL*PRODM
            R_2=R_1/(PRODM*PRODM)
C
            R(-L,M,L)=R_1
            R(-L,-M,L)=R_2
C
            R(L,-M,L)=R_1*CM
            R(M,-L,L)=R_1*CM
            R(L,M,L)=R_2*CM
            R(-M,-L,L)=R_2*CM
C
            R(-M,L,L)=R_1
            R(M,L,L)=R_2
C
            CM2=CM
            DO M2=-L1,-M
              CM2=-CM2
              D0=FLOAT(M2)*FLLM
              CF1=CF(L1,M,-M2)/FLL1
              CF2=FLL1/CF(L,M,-M2)
              IF((M.LT.L1).AND.(-M2.LT.L1)) THEN
                R_A=CF2*((1.-D0)*R(M2,M,L1)-R(M2,M,L-2)*CF1)
                R_B=CF2*((1.+D0)*R(M2,-M,L1)-R(M2,-M,L-2)*CF1)
              ELSE
                R_A=CF2*(1.-D0)*R(M2,M,L1)
                R_B=CF2*(1.+D0)*R(M2,-M,L1)
              ENDIF
C
              R(M2,M,L)=R_A
              R(M2,-M,L)=R_B
C
              R(-M2,-M,L)=R_A*CM2
              R(M,M2,L)=R_A*CM2
              R(-M,M2,L)=R_B*CM2
              R(-M2,M,L)=R_B*CM2
C
              R(-M,-M2,L)=R_A
              R(M,-M2,L)=R_B
C
            ENDDO
          ENDDO
C
          PRODM=PRODM*COEF
          R_1=PRODL*PRODM
          R_2=PRODL/PRODM
          R(-L,L,L)=R_1
          R(L,-L,L)=R_1
          R(L,L,L)=R_2
          R(-L,-L,L)=R_2
C
        ENDDO
      ENDIF
C
      RETURN
C
      END
