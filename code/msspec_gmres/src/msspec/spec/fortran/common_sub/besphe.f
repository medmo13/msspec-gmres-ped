C
C=======================================================================
C
      SUBROUTINE BESPHE(NL,IBES,X1,FL)
C
C  This routine computes the spherical Bessel functions for
C                       a real argument X1.
C
C       IBES=1 : Bessel function
C       IBES=2 : Neumann function
C       IBES=3 : Hankel function of the first kind
C       IBES=4 : Hankel function of the second kind
C       IBES=5 : Modified Bessel function
C       IBES=6 : Modified Neumann function
C       IBES=7 : Modified Hankel function
C
C                                         Last modified : 8 Nov 2006
C
C
      USE DIM_MOD
      USE OUTUNITS_MOD
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 FL(0:2*NL_M),FLNN(0:N_BESS),GL(0:N_BESS),UN,I,C1,C2
      COMPLEX*16 ZERO,CNORM
C
      DOUBLE PRECISION SCALN(0:N_BESS)
C
      REAL X1
C
C
      ECH=37.D0
      COMP=1.D37
      COMM=1.D-37
      X=DBLE(X1)
      NX=INT(X1)
      NREC=5*MAX0(NL-1,NX)
      IF(NREC.GT.N_BESS) GOTO 16
      ITEST=0
      ZERO=(0.D0,0.D0)
      UN=(1.D0,0.D0)
      C1=UN
      I=(0.D0,1.D0)
      C2=I
      DEB=1.D0
      IF((IBES.EQ.3).OR.(IBES.EQ.4)) THEN
        IBES1=1
        IF(IBES.EQ.4) C2=-I
      ELSEIF(IBES.EQ.7) THEN
        IBES1=5
        C2=-UN
      ELSE
        IBES1=IBES
      ENDIF
C
C   Case where the argument is zero
C
      IF(DABS(X).LT.0.000001D0) THEN
        IF((IBES.EQ.1).OR.(IBES.EQ.5)) THEN
          FL(0)=UN
          DO 10 L=1,NL-1
            FL(L)=ZERO
  10      CONTINUE
          ITEST=1
        ELSE
          ITEST=-1
        ENDIF
      ENDIF
      IF(ITEST) 11,12,13
  11  WRITE(IUO1,14)
      STOP
  16  WRITE(IUO1,17) NREC
      STOP
  15  IBES1=IBES1+1
C
C   Initial values
C
  12  A=-1.D0
      B=1.D0
      IF(IBES1.EQ.1) THEN
        FL(0)=UN*DSIN(X)/X
        FLNN(NREC)=ZERO
        SCALN(NREC)=0.D0
        FLNN(NREC-1)=UN*DEB
        SCALN(NREC-1)=0.D0
      ELSEIF(IBES1.EQ.2) THEN
        GL(0)=-UN*DCOS(X)/X
        GL(1)=GL(0)/X -DSIN(X)/X
      ELSEIF(IBES1.EQ.5) THEN
        A=1.D0
        B=-1.D0
        FL(0)=UN*DSINH(X)/X
        FLNN(NREC)=ZERO
        SCALN(NREC)=0.D0
        FLNN(NREC-1)=UN*DEB
        SCALN(NREC-1)=0.D0
      ELSEIF(IBES1.EQ.6) THEN
        A=1.D0
        B=-1.D0
        GL(0)=UN*DCOSH(X)/X
        GL(1)=(DSINH(X)-GL(0))/X
      ENDIF
C
C   Downward reccurence for the spherical Bessel function
C
      IF((IBES1.EQ.1).OR.(IBES1.EQ.5)) THEN
        DO 30 L=NREC-1,1,-1
          ECHEL=0.D0
          SCALN(L-1)=SCALN(L)
          REN=DEXP(SCALN(L)-SCALN(L+1))
          FLNN(L-1)=A*(REN*FLNN(L+1)-B*DFLOAT(2*L+1)*FLNN(L)/X)
          IF(CDABS(FLNN(L-1)).GT.COMP) THEN
            ECHEL=-ECH
          ELSEIF(CDABS(FLNN(L-1)).LT.COMM) THEN
            ECHEL=ECH
          ENDIF
          IF(ECHEL.NE.0.D0 ) SCALN(L-1)=ECHEL+SCALN(L-1)
          FLNN(L-1)=FLNN(L-1)*DEXP(ECHEL)
  30    CONTINUE
        CNORM=FL(0)/FLNN(0)
        DO 40 L=1,NL-1
          FL(L)=CNORM*FLNN(L)*DEXP(SCALN(0)-SCALN(L))
  40    CONTINUE
      ELSE
C
C   Upward recurrence for the spherical Neumann function
C
        DO 20 L=1,NL-1
          IF(IBES.EQ.7) C1=(-UN)**(L+2)
          GL(L+1)=A*GL(L-1)+B*DFLOAT(2*L+1)*GL(L)/X
          IF(IBES1.NE.IBES) THEN
C
C   Calculation of the spherical Hankel function
C
            FL(L+1)=C1*(FL(L+1)+C2*GL(L+1))
          ELSE
            FL(L+1)=GL(L+1)
          ENDIF
  20    CONTINUE
        IF(IBES1.EQ.IBES) THEN
          FL(0)=GL(0)
          FL(1)=GL(1)
        ELSE
          FL(0)=C1*(FL(0)+C2*GL(0))
          FL(1)=C1*(FL(1)+C2*GL(1))
        ENDIF
        IBES1=IBES
      ENDIF
      IF(IBES.NE.IBES1) GOTO 15
C
  13  RETURN
C
  14  FORMAT(/////,3X,'<<<<<<<<<< THE ARGUMENT OF THE BESSEL ','FUNCTION
     &S IS NUL >>>>>>>>>>')
  17  FORMAT(/////,3X,'<<<<<<<<<< THE DIMENSIONNING N_BESS ','IS NOT COR
     &RECT FOR SUBROUTINE BESPHE >>>>>>>>>>',//,15X,'<<<<<<<<<< IT SHOUL
     &D BE AT LEAST : ',I5,' >>>>>>>>>>')
C
      END
