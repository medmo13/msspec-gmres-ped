C
C=======================================================================
C
      SUBROUTINE GAUNT_ST(LMAX_T)
C
C   This subroutine calculates the Gaunt coefficient G(L2,L3|L1)
C    using a downward recursion scheme due to Schulten and Gordon 
C    for the Wigner's 3j symbols. The result is stored as GNT(L3), 
C    making use of the selection rule M3 = M1 - M2.  
C     
C   Ref. :  K. Schulten and R. G. Gordon, J. Math. Phys. 16, 1961 (1975) 
C
C   This is the double precision version where the values are stored
C
C                                          Last modified : 14 May 2009
C
C
      USE DIM_MOD
      USE LOGAMAD_MOD
      USE GAUNT_C_MOD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER LMAX_T
C
      REAL*8 F(0:N_GAUNT),G(0:N_GAUNT),A(0:N_GAUNT),A1(0:N_GAUNT)
      REAL*8 B(0:N_GAUNT)
C
      DATA PI4/12.566370614359D0/
C
      DO L1=0,LMAX_T
        IL1=L1*L1+L1+1
        DO M1=-L1,L1
          IND1=IL1+M1
          LM1=L1+M1
          KM1=L1-M1
          DO L2=0,LMAX_T
            IL2=L2*L2+L2+1
C
            IF(MOD(M1,2).EQ.0) THEN
              COEF=DSQRT(DFLOAT((L1+L1+1)*(L2+L2+1))/PI4)
            ELSE
              COEF=-DSQRT(DFLOAT((L1+L1+1)*(L2+L2+1))/PI4)
            ENDIF
C
            L12=L1+L2
            K12=L1-L2
            L12_1=L12+L12+1
            L12_2=L12*L12
            L12_21=L12*L12+L12+L12+1
            K12_2=K12*K12
C
            F(L12+1)=0.D0
            G(L12+1)=0.D0
            A(L12+1)=0.D0
            A1(L12+1)=0.D0
            A1(L12)=2.D0*DSQRT(DFLOAT(L1*L2*L12_1*L12_2))
            D1=GLD(L2+L2+1,1)-GLD(L12_1+1,1)
            D5=0.5D0*(GLD(L1+L1+1,1)+GLD(L2+L2+1,1)-GLD(L12_1+1,1))
            D6=GLD(L12+1,1)-GLD(L1+1,1)-GLD(L2+1,1)
C
            IF(MOD(K12,2).EQ.0) THEN
              G(L12)=DEXP(D5+D6)
            ELSE
              G(L12)=-DEXP(D5+D6)
            ENDIF
C
            DO M2=-L2,L2
              IND2=IL2+M2
C
              M3=M1-M2
              LM2=L2+M2
              KM2=L2-M2
C
              DO J=1,N_GAUNT
                GNT(J,IND2,IND1)=0.D0
              ENDDO
C
              IF((ABS(M1).GT.L1).OR.(ABS(M2).GT.L2)) GOTO 10
C
              D2=GLD(L1+L1+1,1)-GLD(LM2+1,1)
              D3=GLD(L12+M3+1,1)-GLD(KM2+1,1)
              D4=GLD(L12-M3+1,1)-GLD(LM1+1,1)-GLD(KM1+1,1)
C
              IF(MOD(KM1-KM2,2).EQ.0) THEN
                F(L12)=DSQRT(DEXP(D1+D2+D3+D4))
              ELSE
                F(L12)=-DSQRT(DEXP(D1+D2+D3+D4))
              ENDIF
C
              A(L12)=2.D0*DSQRT(DFLOAT(L1*L2*L12_1*(L12_2-M3*M3)))
              B(L12)=-DFLOAT(L12_1*((L2*L2-L1*L1-K12)*M3+L12*(L12+1)
     1               *(M2+M1)))
C
              IF(ABS(M3).LE.L12) THEN 
                GNT(L12,IND2,IND1)=COEF*F(L12)*G(L12)*
     1                              DSQRT(DFLOAT(L12_1))
              ENDIF
C
              JMIN=MAX0(ABS(K12),ABS(M3))
C
              DO J=L12-1,JMIN,-1
                J1=J+1
                J2=J+2
                JJ=J*J
                A1(J)=DSQRT(DFLOAT(JJ*(JJ-K12_2)*(L12_21-JJ)))
                A(J)=DSQRT(DFLOAT((JJ-K12_2)*(L12_21-JJ)*(JJ-M3*M3)))
                B(J)=-DFLOAT((J+J1)*(L2*(L2+1)*M3-L1*(L1+1)*M3+J*J1*
     1                        (M2+M1)))
                F(J)=-(DFLOAT(J1)*A(J2)*F(J2)+B(J1)*F(J1))/(DFLOAT(J2)*
     1                 A(J1))
                G(J)=-(DFLOAT(J1)*A1(J2)*G(J2))/(DFLOAT(J2)*A1(J1))
C
                IF(ABS(M3).LE.J) THEN 
                  GNT(J,IND2,IND1)=COEF*F(J)*G(J)*DSQRT(DFLOAT(J+J1))
                ENDIF
              ENDDO
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
  10  RETURN
C
      END
 
