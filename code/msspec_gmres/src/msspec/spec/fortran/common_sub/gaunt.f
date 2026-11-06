C
C=======================================================================
C
      SUBROUTINE GAUNT(L2,M2,L1,M1,GNT)
C
C   This subroutine calculates the Gaunt coefficient G(L2,L3|L1)
C    using a downward recursion scheme due to Schulten and Gordon
C    for the Wigner's 3j symbols. The result is stored as GNT(L3),
C    making use of the selection rule M3 = M1 - M2.
C
C   Ref. :  K. Schulten and R. G. Gordon, J. Math. Phys. 16, 1961 (1975)
C
C                                          Last modified :  8 Dec 2008
C
C
      USE DIM_MOD
      USE LOGAMAD_MOD
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL GNT(0:N_GAUNT)
C
      DOUBLE PRECISION F(0:N_GAUNT),G(0:N_GAUNT),A(0:N_GAUNT)
      DOUBLE PRECISION A1(0:N_GAUNT),B(0:N_GAUNT)
C
C
      DATA PI4/12.566370614359D0/
C
      L12=L1+L2
      K12=L1-L2
C
      DO J=1,N_GAUNT
        GNT(J)=0.
      ENDDO
C
      IF((ABS(M1).GT.L1).OR.(ABS(M2).GT.L2)) GOTO 10
C
      M3=M1-M2
      LM1=L1+M1
      LM2=L2+M2
      KM1=L1-M1
      KM2=L2-M2
C
      IF(MOD(M1,2).EQ.0) THEN
        COEF=DSQRT(DFLOAT((2*L1+1)*(2*L2+1))/PI4)
      ELSE
        COEF=-DSQRT(DFLOAT((2*L1+1)*(2*L2+1))/PI4)
      ENDIF
C
      F(L12+1)=0.D0
      G(L12+1)=0.D0
      A(L12+1)=0.D0
      A1(L12+1)=0.D0
      D1=GLD(2*L2+1,1)-GLD(2*L12+2,1)
      D2=GLD(2*L1+1,1)-GLD(LM2+1,1)
      D3=GLD(L12+M3+1,1)-GLD(KM2+1,1)
      D4=GLD(L12-M3+1,1)-GLD(LM1+1,1)-GLD(KM1+1,1)
C
      IF(MOD(KM1-KM2,2).EQ.0) THEN
        F(L12)=DSQRT(DEXP(D1+D2+D3+D4))
      ELSE
        F(L12)=-DSQRT(DEXP(D1+D2+D3+D4))
      ENDIF
C
      D5=0.5D0*(GLD(2*L1+1,1)+GLD(2*L2+1,1)-GLD(2*L12+2,1))
      D6=GLD(L12+1,1)-GLD(L1+1,1)-GLD(L2+1,1)
C
      IF(MOD(K12,2).EQ.0) THEN
        G(L12)=DEXP(D5+D6)
      ELSE
        G(L12)=-DEXP(D5+D6)
      ENDIF
C
      A(L12)=2.D0*DSQRT(DFLOAT(L1*L2*(1+2*L12)*(L12*L12-M3*M3)))
      B(L12)=-DFLOAT((2*L12+1)*((L2*L2-L1*L1-K12)*M3+L12*(L12+1)*(M2+M1)
     &))
      A1(L12)=2.D0*DFLOAT(L12)*DSQRT(DFLOAT(L1*L2*(1+2*L12)))
C
      IF(ABS(M3).LE.L12) THEN
        GNT(L12)=SNGL(COEF*F(L12)*G(L12)*DSQRT(DFLOAT(2*L12+1)))
      ELSE
        GNT(L12)=0.
      ENDIF
C
      JMIN=MAX0(ABS(K12),ABS(M3))
C
      DO J=L12-1,JMIN,-1
        J1=J+1
        J2=J+2
        A(J)=DSQRT(DFLOAT((J*J-K12*K12))*DFLOAT((L12+1)*(L12+1)-J*J)*DFL
     &OAT(J*J-M3*M3))
        B(J)=-DFLOAT((2*J+1)*(L2*(L2+1)*M3-L1*(L1+1)*M3+J*J1*(M2+M1)))
        A1(J)=DFLOAT(J)*DSQRT(DFLOAT((J*J-K12*K12)*((L12+1)*(L12+1)-J*J)
     &))
        F(J)=-(DFLOAT(J1)*A(J2)*F(J2)+B(J1)*F(J1))/(DFLOAT(J2)*A(J1))
        G(J)=-(DFLOAT(J1)*A1(J2)*G(J2))/(DFLOAT(J2)*A1(J1))
        GND=COEF*F(J)*G(J)*DSQRT(DFLOAT(2*J+1))
C
        IF(ABS(M3).LE.J) THEN
          GNT(J)=SNGL(GND)
        ELSE
          GNT(J)=0.
        ENDIF
C
      ENDDO
C
  10  RETURN
C
      END
