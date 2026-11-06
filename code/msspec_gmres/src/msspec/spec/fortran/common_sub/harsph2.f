C
C=======================================================================
C
      SUBROUTINE HARSPH2(NL,THETA,PHI,YLM,NC)
C
C  This routine computes the complex spherical harmonics using Condon and
C                  Shortley phase convention. This version for m=0 only
C
      USE DIM_MOD
C
      USE EXPFAC2_MOD
      USE FACTSQ_MOD
C
      COMPLEX YLM(0:NL,-NL:NL),COEF,YMM,YMMP,C
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
      RETURN
C
      END
