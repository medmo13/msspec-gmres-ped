C
C=======================================================================
C
      SUBROUTINE N_J(J1,MJ1,J2,MJ2,MJ6,NJ,I_INT,N_IN)
C
C   This subroutine calculates Wigner's 3j and 6j coefficients
C    using a downward recursion scheme due to Schulten and Gordon.
C    The 3j are defined as (J1 J2 J) where in fact L1=MJ1, etc are
C                          (L1 L2 L)
C    azimuthal quantum numbers, and the 6j as {J1 J2 J} where now
C                                             {L1 L2 L}
C    J1, L1, etc are the same kind of orbital quantum numbers.
C    The result is stored as NJ(J).
C
C    The parameter N allows to choose between 3j and 6j calculation, and
C    Clebsch-Gordan. It can take the values :
C
C                    N = 2 ----> Clebsch-Gordan
C                    N = 3 ----> Wigner's 3j
C                    N = 6 ----> Wigner's 6j
C
C    The Clebsch-Gordan coefficients are related to Wigner's 3j through :
C
C       CG(J1,M1,J2,M2|J,MJ) = ( J1 J2  J  )*sqrt(2*J+1)*(-1)**(J1-J2+MJ)
C                              ( M1 M2 -MJ )
C    I_INT is a flag that returns 1 if the index J of the nj symbol
C      is integer and 0 if it is a half integer.
C
C   Note : For 3j, MJ6 is ignored while for 6j, we have :
C
C                J1=J1   MJ1=L1   J2=J2   MJ2=L2   MJ6=L
C
C   Ref. : K. Schulten and R. G. Gordon, J. Math. Phys. 16, 1961 (1975)
C
C                      Last modified :  8 Dec 2008 ----> D. Sebilleau
C
C
      USE DIM_MOD
      USE LOGAMAD_MOD
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*4 NJ(0:N_GAUNT)
C
      REAL*8 J1,J2,J,MJ1,MJ2,MJ,JP1,JP2
      REAL*8 F(0:N_GAUNT),A(0:N_GAUNT),B(0:N_GAUNT)
      REAL*8 JL12,JK12,MJ6,SIG
      REAL*8 JJ1,JJ2,JL1,JL2,JL3,JJ12,JJ_MIN,JJ_MAX
C
C
      DATA SMALL /0.0001D0/
C
      IS=0
C
      IF(N_IN.EQ.2) THEN
        N_OU=3
        I_CG=1
      ELSE
        N_OU=N_IN
        I_CG=0
      ENDIF
C
      IF(N_OU.EQ.3) THEN
C
C------------------------------  3j case  ---------------------------------
C
C
C  Test to check if J1 and J2 are integer or semi-integer
C
C     Integer      : IG=1
C     Half-integer : IG=2
C
C  Each angular momentum J is represented by the integer index L and
C    the corresponding MJ by M
C
        L1=INT(J1+SMALL)
        L2=INT(J2+SMALL)
        M1=INT(MJ1+SIGN(SMALL,MJ1))
        M2=INT(MJ2+SIGN(SMALL,MJ2))
        DIF1=J1-DFLOAT(L1)
        DIF2=J2-DFLOAT(L2)
C
C  IGx is a flag telling the code which case of Gamma function to use :
C
C     IGx = 1 : integer case
C     IGx = 2 : half-integer case
C
        IF(ABS(DIF1).LT.SMALL) THEN
          IG1=1
        ELSE
          IG1=2
        ENDIF
        IF(ABS(DIF2).LT.SMALL) THEN
          IG2=1
        ELSE
          IG2=2
        ENDIF
        IF(IG1.EQ.IG2) THEN
          IGG=1
          IF(IG1.EQ.2) IS=1
        ELSE
          IGG=2
        ENDIF
C
C  Here, we assume that (J1,J2) are both either integer or half-integer
C  If J is integer, the corresponding index is L = j (for loops or storage)
C    while if J is an half-integer, this index is L=  j - 1/2 = int(j)
C
C  Integer indices are used for loops and for storage while true values
C    are used for the initial values. When J1 and J2 are both half-integers,
C    the values of J are integer and L should be increased by 1
C
        JL12=J1+J2
        JK12=J1-J2
C
        L12=INT(JL12 + SIGN(SMALL,JL12))
        K12=INT(JK12 + SIGN(SMALL,JK12))
C
        LM1=INT(J1+MJ1 + SIGN(SMALL,J1+MJ1))
        LM2=INT(J2+MJ2 + SIGN(SMALL,J2+MJ2))
        KM1=INT(J1-MJ1 + SIGN(SMALL,J1-MJ1))
        KM2=INT(J2-MJ2 + SIGN(SMALL,J2-MJ2))
C
        MJ=-MJ1-MJ2
C
        M=INT(MJ+SIGN(SMALL,MJ))
        L12M=INT(JL12+MJ+SIGN(SMALL,JL12+MJ))
        K12M=INT(JL12-MJ+SIGN(SMALL,JL12-MJ))
        L1_2=INT(J1+J1+SIGN(SMALL,J1))
        L2_2=INT(J2+J2+SIGN(SMALL,J2))
        L12_2=INT(JL12+JL12+SIGN(SMALL,JL12))
C
        IF(IG(JL12).EQ.1) THEN
          I_INT=1
        ELSE
          I_INT=0
        ENDIF
C
C  Initialisation of the 3j symbol NJ(J) = (J1  J2   J)
C                                          (MJ1 MJ2 MJ)
C
        DO L=0,L12
          NJ(L)=0.
        ENDDO
C
        IF((ABS(MJ1).GT.J1).OR.(ABS(MJ2).GT.J2)) GOTO 10
C
C  Initial values (J1+J2+1) and (J1+J2) for J to be used in the downward
C     recursion scheme. This scheme writes as
C
C     J A(J+1) NJ(J+1) + B(J) NJ(J) + (J+1) A(J) NJ(J-1) = 0
C
        F(L12+1)=0.D0
        A(L12+1)=0.D0
        D1=GLD(L2_2+1,1)-GLD(L12_2+2,1)
        D2=GLD(L1_2+1,1)-GLD(LM2+1,1)
        D3=GLD(L12M+1,1)-GLD(KM2+1,1)
        D4=GLD(K12M+1,1)-GLD(LM1+1,1)-GLD(KM1+1,1)
C
        N12=INT(JK12-MJ + SIGN(SMALL,JK12-MJ))
C
        IF(I_CG.EQ.1) THEN
          IF(MOD(N12,2).EQ.0) THEN
            SIG=1.D0
          ELSE
            SIG=-1.D0
          ENDIF
        ENDIF
C
        IF(MOD(N12,2).EQ.0) THEN
          F(L12)=DSQRT(DEXP(D1+D2+D3+D4))
        ELSE
          F(L12)=-DSQRT(DEXP(D1+D2+D3+D4))
       ENDIF
C
        A(L12)=2.D0*DSQRT(J1*J2*(1.D0+2.D0*JL12)*(JL12*JL12-MJ*MJ))
        B(L12)=-(2.D0*JL12+1.D0)*((J1*J1-J2*J2+JK12)*MJ-JL12*(JL12+1.D0)
     &*(MJ2-MJ1))
C
        IF(ABS(M).LE.L12) THEN
          IF(I_CG.EQ.0) THEN
            NJ(L12)=SNGL(F(L12))
          ELSE
            NJ(L12)=SNGL(F(L12)*SIG*DSQRT(JL12+JL12+1.D0))
          ENDIF
        ELSE
          NJ(L12)=0.
        ENDIF
C
        LMIN=MAX0(ABS(K12),ABS(M))
C
C  Downward recursion for NJ(J)
C
        DO L=L12-1,LMIN,-1
          LP1=L+1
          LP2=L+2
C
C  Value of the angular momentum J corresponding to the loop index L
C
          IF(IGG.EQ.1) THEN
            J=DFLOAT(L)
            JP1=DFLOAT(LP1)
            JP2=DFLOAT(LP2)
          ELSE
            J=DFLOAT(L) + 0.5D0
            JP1=DFLOAT(LP1) + 0.5D0
            JP2=DFLOAT(LP2) + 0.5D0
          ENDIF
C
          A(L)=DSQRT((J*J-JK12*JK12)*((JL12+1.D0)*(JL12+1.D0)-J*J)*(J*J-
     &MJ*MJ))
          B(L)=-(2.D0*J+1.D0)*(J1*(J1+1.D0)*MJ-J2*(J2+1.D0)*MJ-J*JP1*(MJ
     &2-MJ1))
          F(L)=-(JP1*A(LP2)*F(LP2)+B(LP1)*F(LP1))/(JP2*A(LP1))
C
          IF(ABS(MJ).LE.J) THEN
            IF(I_CG.EQ.0) THEN
              NJ(L)=SNGL(F(L))
            ELSE
              NJ(L)=SNGL(F(L)*SIG*DSQRT(J+J+1.D0))
            ENDIF
          ELSE
            NJ(L)=0.
          ENDIF
C
        ENDDO
C
  10    CONTINUE
C
      ELSEIF(N_OU.EQ.6) THEN
C
C------------------------------  6j case  ---------------------------------
C
C  Change of notation for greater readability ---> NJ(JJ)
C
C    True angular momentum value : begins with a J (JJn,JLn)
C    Corresponding integer storage and loop index : begins by L (LJn,LLn)
C
        JJ1=J1
        JJ2=J2
        JL1=MJ1
        JL2=MJ2
        JL3=MJ6
C
        LJ1=INT(JJ1+SIGN(SMALL,JJ1))
        LJ2=INT(JJ2+SIGN(SMALL,JJ2))
        LL1=INT(JL1+SIGN(SMALL,JL1))
        LL2=INT(JL2+SIGN(SMALL,JL2))
        LL3=INT(JL3+SIGN(SMALL,JL3))
C
        JJ12=JJ1-JJ2
        JL12=JL1-JL2
C
        LJ12=INT(JJ12+SIGN(SMALL,JJ12))
        LL12=INT(JL12+SIGN(SMALL,JL12))
C
        JJ_MIN=MAX(ABS(LJ12),ABS(LL12))
        JJ_MAX=MIN(JJ1+JJ2,JL1+JL2)
        LJJ_MIN=INT(JJ_MIN+SIGN(SMALL,JJ_MIN))
        LJJ_MAX=INT(JJ_MAX+SIGN(SMALL,JJ_MAX))
C
C  Initialisation of the 6j symbol NJ(J) = {J1  J2  J }
C                                          {L1  L2  L3}
C
        DO L=0,LJJ_MAX
          NJ(L)=0.
        ENDDO
C
C  Initial values (J1+J2+1) and (J1+J2) for J to be used in the downward
C     recursion scheme. This scheme writes as
C
C     J A(J+1) NJ(J+1) + B(J) NJ(J) + (J+1) A(J) NJ(J-1) = 0
C
C  There are two possible initial values as max(|J1-J2|,|L1-L2|) <= J <=
C    min(J1+J2,L1+L2) :
C
C    {J1  J2  L1+L2}  and  {J1  J2  J1+J2}  =  {L1 L2 J1+J2}
C    {L1  L2    L3 }       {L1  L2    L3 }     {J1 J2   L3 }
C
C    They can be calculated from equation (6.3.1) of Edmonds page 97
C
        F(LJJ_MAX+1)=0.D0
        A(LJJ_MAX+1)=0.D0
C
        IF(ABS(JJ_MAX-JL1-JL2).LT.SMALL) THEN
          F(LJJ_MAX)=SIXJ_IN(JJ1,JJ2,JL1,JL2,JL3)
        ELSE
          F(LJJ_MAX)=SIXJ_IN(JL1,JL2,JJ1,JJ2,JL3)
        ENDIF
        NJ(LJJ_MAX)=SNGL(F(LJJ_MAX))
C
        A(LJJ_MAX)=SQRT((JJ_MAX*JJ_MAX-(JJ1-JJ2)*(JJ1-JJ2))*((JJ1+JJ2+1.
     &D0)*(JJ1+JJ2+1.D0)-JJ_MAX*JJ_MAX)*(JJ_MAX*JJ_MAX-(JL1-JL2)*(JL1-JL
     &2))*((JL1+JL2+1.D0)*(JL1+JL2+1.D0)-JJ_MAX*JJ_MAX))
        B(LJJ_MAX)=(JJ_MAX+JJ_MAX+1.D0)*(JJ_MAX*(JJ_MAX+1.D0)*(-JJ_MAX*(
     &JJ_MAX+1.D0)+JJ1*(JJ1+1.D0)+JJ2*(JJ2+1.D0))+JL1*(JL1+1.D0)*(JJ_MAX
     &*(JJ_MAX+1.D0)+JJ1*(JJ1+1.D0)-JJ2*(JJ2+1.D0))+JL2*(JL2+1.D0)*(JJ_M
     &AX*(JJ_MAX+1.D0)-JJ1*(JJ1+1.D0)+JJ2*(JJ2+1.D0))-(JJ_MAX+JJ_MAX)*(J
     &J_MAX+1.D0)*JL3*(JL3+1.D0))
C
      IF(IG(JJ_MAX).EQ.1) THEN
        I_INT=1
      ELSE
        I_INT=0
      ENDIF
C
C  Downward recurrence relation
C
        DO L=LJJ_MAX-1,LJJ_MIN,-1
          LP1=L+1
          LP2=L+2
C
C  Value of the angular momentum J corresponding to the loop index L
C
          IF(IG(JJ_MAX).EQ.1) THEN
            J=DFLOAT(L)
            JP1=DFLOAT(LP1)
            JP2=DFLOAT(LP2)
          ELSE
            J=DFLOAT(L) + 0.5D0
            JP1=DFLOAT(LP1) + 0.5D0
            JP2=DFLOAT(LP2) + 0.5D0
          ENDIF
C
          A(L)=SQRT((J*J-(JJ1-JJ2)*(JJ1-JJ2))*((JJ1+JJ2+1.D0)*(JJ1+JJ2+1
     &.D0)-J*J)*(J*J-(JL1-JL2)*(JL1-JL2))*((JL1+JL2+1.D0)*(JL1+JL2+1.D0)
     &-J*J))
          B(L)=(J+J+1)*(J*JP1*(-J*JP1+JJ1*(JJ1+1.D0)+JJ2*(JJ2+1.D0))+JL1
     &*(JL1+1.D0)*(J*JP1+JJ1*(JJ1+1.D0)-JJ2*(JJ2+1.D0))+JL2*(JL2+1.D0)*(
     &J*JP1-JJ1*(JJ1+1.D0)+JJ2*(JJ2+1.D0))-(J+J)*JP1*JL3*(JL3+1.D0))
C
          F(L)=-(JP1*A(LP2)*F(LP2)+B(LP1)*F(LP1))/(JP2*A(LP1))
          NJ(L)=SNGL(F(L))
C
        ENDDO
C
      ENDIF
C
      RETURN
C
      END
