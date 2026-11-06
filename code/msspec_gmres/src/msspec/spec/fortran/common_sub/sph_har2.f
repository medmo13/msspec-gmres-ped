C
C=======================================================================
C
      SUBROUTINE SPH_HAR2(NL,X,CF,YLM,NC)
C
C  This routine computes the complex spherical harmonics using Condon and
C                  Shortley phase convention.
C
C  If the angular direction R=(THETAR,PHIR)  is given in cartesian
C      coordinates by (XR,YR,ZR), the arguments of the subroutine are :
C
C                    X  = ZR         = cos(THETAR)
C                    CF = XR + i YR  = sin(THETAR)*exp(i PHIR)
C
C          NL is the dimensioning of the YLM array and NC is
C                   the maximum l value to be computed.
C
C  This is the double precision version of sph_har.f
C
C
      USE DIM_MOD
C
      USE DEXPFAC2_MOD
      USE DFACTSQ_MOD
C
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 YLM(0:NL,-NL:NL),COEF,YMM,YMMP,C,CF
C
      DATA SQ4PI_INV,SQR3_INV /0.282094791774D0,0.488602511903D0/
C
C
      YLM(0,0)=DCMPLX(SQ4PI_INV)
      YLM(1,0)=X*SQR3_INV
      DO L=2,NC
        Y=1.D0/DFLOAT(L)
        YLM(L,0)=X*DSQRT(4.D0-Y*Y)*YLM(L-1,0) - (1.D0-Y)*DSQRT(1.D0+2.D0
     &/(DFLOAT(L)-1.5D0))*YLM(L-2,0)
      ENDDO
C
      C2=-1.D0
      C=-0.5D0*CF
C
      C1=1.D0
      COEF=(1.D0,0.D0)
      DO M=1,NC
        C1=C1*C2
        COEF=COEF*C
        YMM=SQ4PI_INV*COEF*DFSQ(M)
        YLM(M,M)=YMM
        YLM(M,-M)=C1*DCONJG(YMM)
        YMMP=X*DSQRT(DFLOAT(M+M+3))*YMM
        YLM(M+1,M)=YMMP
        YLM(M+1,-M)=C1*DCONJG(YMMP)
        IF(M.LT.NC-1) THEN
          DO L=M+2,NC
            YLM(L,M)=(X*DFLOAT(L+L-1)*DEXPF2(L-1,M)*YLM(L-1,M) - DFLOAT(
     &L+M-1)*DEXPF2(L-2,M)*YLM(L-2,M))/(DEXPF2(L,M)*DFLOAT(L-M))
            YLM(L,-M)=C1*DCONJG(YLM(L,M))
          ENDDO
        ENDIF
      ENDDO
C
      RETURN
C
      END
