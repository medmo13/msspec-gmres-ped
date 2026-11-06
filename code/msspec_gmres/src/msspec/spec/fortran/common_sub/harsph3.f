C
C=======================================================================
C
      SUBROUTINE HARSPH3(NL,THETA,PHI,YLM2,NC)
C
C  This routine computes the complex spherical harmonics using Condon and
C                  Shortley phase convention.
C
      USE DIM_MOD
C
      USE EXPFAC2_MOD
      USE FACTSQ_MOD
C
      COMPLEX YLM(0:NL,-NL:NL),COEF,YMM,YMMP,C
      COMPLEX YLM2(LINMAX)
C
      DATA SQ4PI_INV,SQR3_INV /0.282095,0.488602/
      DATA PI,SMALL /3.141593,0.0001/
C
      X=COS(THETA)
      IF(ABS(X).LT.SMALL) X=0.0
      IF(ABS(X+1.).LT.SMALL) X=-1.0
      IF(ABS(X-1.).LT.SMALL) X=1.0
C
      YLM(0,0)=CMPLX(SQ4PI_INV)
      YLM(1,0)=X*SQR3_INV
      DO L=2,NC
        Y=1./FLOAT(L)
        YLM(L,0)=X*SQRT(4.-Y*Y)*YLM(L-1,0) - (1.-Y)*SQRT(1.+2./(FLOAT(L)
     &-1.5))*YLM(L-2,0)
      ENDDO
C
      C2=-1.
      IF((THETA.GE.0.).AND.(THETA.LE.PI)) THEN
        C=-0.5*SQRT(1.-X*X)*EXP((0.,1.)*PHI)
      ELSE
        C=0.5*SQRT(1.-X*X)*EXP((0.,1.)*PHI)
      ENDIF
C
      C1=1.
      COEF=(1.,0.)
      DO M=1,NC
        C1=C1*C2
        COEF=COEF*C
        YMM=SQ4PI_INV*COEF*FSQ(M)
        YLM(M,M)=YMM
        YLM(M,-M)=C1*CONJG(YMM)
        YMMP=X*SQRT(FLOAT(M+M+3))*YMM
        YLM(M+1,M)=YMMP
        YLM(M+1,-M)=C1*CONJG(YMMP)
        IF(M.LT.NC-1) THEN
          DO L=M+2,NC
            YLM(L,M)=(X*(L+L-1)*EXPF2(L-1,M)*YLM(L-1,M) - (L+M-1)*EXPF2(
     &L-2,M)*YLM(L-2,M))/(EXPF2(L,M)*(L-M))
            YLM(L,-M)=C1*CONJG(YLM(L,M))
          ENDDO
        ENDIF
      ENDDO
C
      DO L=0,NC
        IL=L*L+L+1
        DO M=-L,L
          IND=IL+M
          YLM2(IND)=YLM(L,M)
        ENDDO
      ENDDO
C
      RETURN
C
      END
