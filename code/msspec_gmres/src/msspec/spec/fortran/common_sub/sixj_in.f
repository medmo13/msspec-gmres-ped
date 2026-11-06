C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SIXJ_IN(J1,J2,L1,L2,L3)
C
C  This function calculates the initial value {J1 J2 L1+L2}
C                                             {L1 L2   L3 }
C
C  A 6j symbol {J1 J2 J3} is non zero only if
C              {J4 J5 J6}
C
C   (J1,J2,J3),(J4,J5,J3),(J2,J4,J6) and (J1,J5,J6) satisfy the triangular inequality :
C
C       (a,b,c) non zero if |a-b| <= c <= (a+b) . This means also that (a+b) and c must
C       have the same nature (integer or half-integer).
C
C   (J1,J2,J3) and (J4,J5,J3) are taken care of by the bounds of J3, JJ_MIN and JJ_MAX,
C       as chosen in the N_J routine. Here we check the two last ones.
C
C                                        Last modified :  8 Dec 2008
C
C
      USE DIM_MOD
      USE LOGAMAD_MOD
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 J1,J2,L1,L2,L3
C
C
      DATA SMALL /0.0001/
C
      IZERO=0
C
C  Check for unphysical values of L3
C
      IF(DABS(J2-L1).GT.L3) IZERO=1
      IF(J2+L1.LT.L3) IZERO=1
      IF(IG(J2+L1).NE.IG(L3)) IZERO=1
      IF(DABS(J1-L2).GT.L3) IZERO=1
      IF(J1+L2.LT.L3) IZERO=1
      IF(IG(J1+L2).NE.IG(L3)) IZERO=1
C
      IF(IZERO.EQ.1) THEN
        SIXJ_IN=0.D0
      ELSE
C
C  Storage indices of the angular momenta.
C
        LJ1=INT(J1+SIGN(SMALL,J1))
        LJ2=INT(J2+SIGN(SMALL,J2))
        LL1=INT(L1+SIGN(SMALL,L1))
        LL2=INT(L2+SIGN(SMALL,L2))
        LL3=INT(L3+SIGN(SMALL,L3))
        LL1_2=INT(L1+L1+SIGN(SMALL,L1))
        LL2_2=INT(L2+L2+SIGN(SMALL,L2))
C
        MSIGN=INT(J1+J2+L1+L2+SIGN(SMALL,J1+J2+L1+L2))
        IF(MOD(MSIGN,2).EQ.0) THEN
          SIGNE=1.D0
        ELSE
          SIGNE=-1.D0
        ENDIF
C
        D1=GLD(LL1_2+1,1) + GLD(LL2_2+1,1) - GLD(LL1_2+LL2_2+2,1)
        D2=GLD(INT(J1+J2+L1+L2)+2,IG(J1+J2+L1+L2)) - GLD(INT(J1+J2-L1-L2
     &)+1,IG(J1+J2-L1-L2))
        D3=GLD(INT(J1-J2+L1+L2)+1,IG(J1-J2+L1+L2)) - GLD(INT(J1+L2-L3)+1
     &,IG(J1+L2-L3))
        D4=GLD(INT(J2-J1+L1+L2)+1,IG(J2-J1+L1+L2)) -GLD(INT(-J1+L2+L3)+1
     &,IG(-J1+L2+L3))
        D5=GLD(INT(J1-L2+L3)+1,IG(J1-L2+L3)) - GLD(INT(J1+L2+L3)+2,IG(J1
     &+L2+L3))
        D6=GLD(INT(J2+L3-L1)+1,IG(J2+L3-L1)) - GLD(INT(J2-L3+L1)+1,IG(J2
     &-L3+L1))
        D7=GLD(INT(L1+L3-J2)+1,IG(L1+L3-J2)) +GLD(INT(L1+J2+L3)+2,IG(L1+
     &J2+L3))
C
        SIXJ_IN=SIGNE*DSQRT(DEXP(D1+D2+D3+D4+D5+D6-D7))
C
      ENDIF
C
      END
