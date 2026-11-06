C
C=======================================================================
C
      SUBROUTINE PATHOP(JPOS,JORDP,JE,I_CP,RHO01,PHI01,RHOIJ,THIJ,PHIIJ,
     &FREF,IJ,D,TAU)
C
C  This subroutine calculates the contribution of a given path to
C    the scattering path operator TAU.
C
C                                           Last modified :  3 Aug 2007
C
      USE DIM_MOD
      USE PATH_SAMPLING_MOD
C
      USE APPROX_MOD
      USE C_RENORM_MOD
      USE EXPFAC_MOD
      USE EXTREM_MOD
      USE INIT_L_MOD
      USE INIT_J_MOD
      USE LBD_MOD
      USE LINLBD_MOD
      USE OUTUNITS_MOD
      USE PATH_MOD
      USE PRINTP_MOD
      USE RA_MOD
      USE RENORM_MOD
      USE ROT_MOD
      USE SCATMAT_MOD , F => F21
      USE TESTS_MOD
      USE TLDW_MOD
      USE TRANS_MOD
      USE VARIA_MOD
C
      LOGICAL MC_SKIP
      DOUBLE PRECISION MC_WEIGHT
      INTEGER JPOS(NDIF_M,3),AMU1
C
C
      REAL RLMIJ(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1)
C
      COMPLEX TAU(LINMAX,LINFMAX,NATCLU_M)
      COMPLEX H(NLAMBDA_M,NLAMBDA_M)
      COMPLEX G(NLAMBDA_M,NLAMBDA_M)
      COMPLEX HLM01(0:NO_ST_M,0:NL_M-1),HLMIJ(0:NO_ST_M,0:NL_M-1)
      COMPLEX SUM_NUJ_0,SUM_MUJ_0,SUM_NU1_0
      COMPLEX SUM_NUJ_1,SUM_MUJ_1,SUM_NU1_1
      COMPLEX SUM_NU1_2,SUM_NU1_3
      COMPLEX RHO01,RHOIJ
      COMPLEX RLMF_0,RLMF_1
      COMPLEX CF,CJ,OVK
      COMPLEX EXP_J,EXP_F,SUM_1
      COMPLEX TL_J
      COMPLEX COEF,ONEC,ZEROC
C
C
C
      DATA PI,XCOMP /3.141593,1.E-10/
C
      CALL SAMPLE_PATH(JORDP,MC_WEIGHT,MC_SKIP)
      IF(MC_SKIP) RETURN
C
      ZEROC=(0.,0.)
      ONEC=(1.,0.)
C
      OVK=(1.,0.)/VK(JE)
      IF(NPATHP.GT.0) THEN
        FM1=FMIN(JORDP)
        XMAX=0.
      ENDIF
      EXP_J=CEXP((0.,-1.)*(PHIIJ-PI))
      EXP_F=CEXP((0.,1.)*PHI01)
      JTYP=JPOS(JORDP,1)
      ITYP=JPOS(1,1)
      JATL=JPOS(JORDP,3)
      IF(I_CP.EQ.0) THEN
        LMJ=LMAX(JTYP,JE)
      ELSE
        LMJ=LF2
      ENDIF
      IF(NO.EQ.8) THEN
        NN2=LMJ+1
      ELSE
        NN2=NO
      ENDIF
      IF(NO.GT.LF2) THEN
        NN=LF2
      ELSE
        NN=NO
      ENDIF
C
C  NO is atom-dependent and is decreased with the rank of the scatterer
C    in the path when I_NO > 0 (except for the first scatterer ITYP for
C    which there is no such decrease)
C
      NO1=N_RA(ITYP)
      IF(I_NO.EQ.0) THEN
        IF(IJ.EQ.1) THEN
          NOJ=N_RA(JTYP)
        ELSE
          NOJ=0
        ENDIF
      ELSE
        IF(IJ.EQ.1) THEN
          NOJ= MAX(N_RA(JTYP)-(JORDP-1)/I_NO,0)
        ELSE
          NOJ=0
        ENDIF
      ENDIF
      NUMX=NO1/2
      NUMAXJ=NOJ/2
C
C   Calculation of the attenuation coefficients along the path
C
      COEF=CEX(1)*OVK
      COEF=COEF*MC_WEIGHT
      DO JSC=2,JORDP
        COEF=COEF*CEXDW(JSC)
      ENDDO
C
C  Renormalization of the path
C
      IF(I_REN.GE.1) THEN
         COEF=COEF*C_REN(JORDP)
C         write(354,*) JORDP,C_REN(JORDP)
      ENDIF
C
C   Call of the subroutines used for the R-A termination matrix
C      This termination matrix is now merged into PATHOP
C
      CALL DJMN2(-THIJ,RLMIJ,LMJ,1)
      CALL POLHAN(ISPHER,NN,LF2,RHO01,HLM01)
      CALL POLHAN(ISPHER,NN2,LMJ,RHOIJ,HLMIJ)
C
      LBD1M1=LBDM(1,1)
      LBD1M2=LBDM(2,1)
C
C   Calculation of the L-independent part of TAU, called H
C
      IF(JORDP.GE.3) THEN
        DO JPAT=2,JORDP-1
          LBD2M=LBDM(1,JPAT)
          LBD3M=LBDM(2,JPAT)
          DO LAMBDA1=1,LBD1M1
            DO LAMBDA3=1,LBD3M
              SUM_1=ZEROC
              DO LAMBDA2=1,LBD2M
                IF(JPAT.GT.2) THEN
                  SUM_1=SUM_1+H(LAMBDA2,LAMBDA1)*F(1,LAMBDA3,LAMBDA2,JPA
     &T)
                ELSE
                  SUM_1=SUM_1+F(1,LAMBDA2,LAMBDA1,1)*F(1,LAMBDA3,LAMBDA2
     &,2)
                ENDIF
              ENDDO
              G(LAMBDA3,LAMBDA1)=SUM_1
            ENDDO
          ENDDO
          DO LAMBDA1=1,LBD1M1
            DO LAMBDA2=1,LBD3M
              H(LAMBDA2,LAMBDA1)=G(LAMBDA2,LAMBDA1)
            ENDDO
          ENDDO
        ENDDO
      ELSEIF(JORDP.EQ.2) THEN
        DO LAMBDA1=1,LBD1M1
          DO LAMBDA2=1,LBD1M2
            H(LAMBDA2,LAMBDA1)=F(1,LAMBDA2,LAMBDA1,1)
          ENDDO
        ENDDO
      ELSEIF(JORDP.EQ.1) THEN
        DO LAMBDA1=1,LBD1M1
          DO LAMBDA2=1,LBD1M1
            H(LAMBDA2,LAMBDA1)=ONEC
          ENDDO
        ENDDO
      ENDIF
C
C   Calculation of the path operator TAU
C
      DO LF=LF1,LF2,ISTEP_LF
        ILF=LF*LF+LF+1
C
        NU1MAX1=MIN(LF,NUMX)
C
C  Case MF = 0
C
        DO LJ=0,LMJ
          ILJ=LJ*LJ+LJ+1
          NUJMAX=MIN(LJ,NUMAXJ)
          IF(JORDP.EQ.1) THEN
            NU1MAX=MIN(NU1MAX1,LJ)
          ELSE
            NU1MAX=NU1MAX1
          ENDIF
C
          IF(ISPEED.EQ.1) THEN
            TL_J=COEF*TL(LF,1,1,JE)*TL(LJ,1,JTYP,JE)
          ELSE
            TL_J=COEF*TLT(LF,1,1,JE)*TLT(LJ,1,JTYP,JE)
          ENDIF
C
C    Case MJ = 0
C
          SUM_NU1_0=ZEROC
C
          DO NU1=0,NU1MAX
            IF(JORDP.GT.1) THEN
              MU1MAX=MIN(LF-NU1,NO1-NU1-NU1)
            ELSE
              MU1MAX=MIN(LF-NU1,NO1-NU1-NU1,LJ)
            ENDIF
C
            DO MU1=-MU1MAX,MU1MAX
              LAMBDA1=LBD(MU1,NU1)
              AMU1=ABS(MU1)
C
              RLMF_0=HLM01(AMU1+NU1,LF)*RLM01(MU1,0,LF)
C
              SUM_NUJ_0=ZEROC
C
              IF(JORDP.GT.1) THEN
                DO NUJ=0,NUJMAX
                  MUJMAX=MIN(LJ,NOJ-NUJ-NUJ)
C
                  SUM_MUJ_0=ZEROC
C
                  DO MUJ=-MUJMAX,MUJMAX
C
                    LAMBDAJ=LBD(MUJ,NUJ)
C
                    SUM_MUJ_0=SUM_MUJ_0+H(LAMBDAJ,LAMBDA1)*RLMIJ(MUJ,0,L
     &J)
                  ENDDO
                  SUM_NUJ_0=SUM_NUJ_0+SUM_MUJ_0*HLMIJ(NUJ,LJ)
C
                ENDDO
              ELSE
                SUM_NUJ_0=HLMIJ(NU1,LJ)*RLMIJ(MU1,0,LJ)
              ENDIF
C
              SUM_NU1_0=SUM_NU1_0+RLMF_0*SUM_NUJ_0
C
            ENDDO
C
          ENDDO
C
          TAU(ILJ,ILF,JATL)=TAU(ILJ,ILF,JATL)+TL_J*SUM_NU1_0
C
          IF(NPATHP.EQ.0) GOTO 35
C
          FM2=FMAX(JORDP)
          XINT=CABS(TL_J*SUM_NU1_0)
          XMAX=AMAX1(XINT,XMAX)
          FMAX(JORDP)=AMAX1(FM2,XINT)
          IF((FMAX(JORDP)-FM2).GT.XCOMP) NPMA(JORDP)=NPATH(JORDP)
          IF((IREF.EQ.1).AND.(JORDP.EQ.NCUT)) THEN
              FREF=FMAX(JORDP)
          ENDIF
 35       CONTINUE
C
C    Case MJ > 0
C
          CJ=ONEC
          DO MJ=1,LJ
            INDJ=ILJ+MJ
            INDJP=ILJ-MJ
            CJ=CJ*EXP_J
C
            SUM_NU1_0=ZEROC
            SUM_NU1_1=ZEROC
C
            DO NU1=0,NU1MAX
              IF(JORDP.GT.1) THEN
                MU1MAX=MIN(LF-NU1,NO1-NU1-NU1)
              ELSE
                MU1MAX=MIN(LF-NU1,NO1-NU1-NU1,LJ)
              ENDIF
C
              DO MU1=-MU1MAX,MU1MAX
                LAMBDA1=LBD(MU1,NU1)
                AMU1=ABS(MU1)
C
                RLMF_0=HLM01(AMU1+NU1,LF)*RLM01(MU1,0,LF)
C
                SUM_NUJ_0=ZEROC
                SUM_NUJ_1=ZEROC
C
                IF(JORDP.GT.1) THEN
                  DO NUJ=0,NUJMAX
                    MUJMAX=MIN(LJ,NOJ-NUJ-NUJ)
C
                    SUM_MUJ_0=ZEROC
                    SUM_MUJ_1=ZEROC
C
                    DO MUJ=-MUJMAX,MUJMAX
C
                      LAMBDAJ=LBD(MUJ,NUJ)
C
                      SUM_MUJ_1=SUM_MUJ_1+H(LAMBDAJ,LAMBDA1)*RLMIJ(MUJ,-
     &MJ,LJ)
                      SUM_MUJ_0=SUM_MUJ_0+H(LAMBDAJ,LAMBDA1)*RLMIJ(MUJ,M
     &J,LJ)
C
                    ENDDO
C
                    SUM_NUJ_0=SUM_NUJ_0+SUM_MUJ_0*HLMIJ(NUJ,LJ)
                    SUM_NUJ_1=SUM_NUJ_1+SUM_MUJ_1*HLMIJ(NUJ,LJ)
C
                  ENDDO
                ELSE
                  SUM_NUJ_1=HLMIJ(NU1,LJ)*RLMIJ(MU1,-MJ,LJ)
                  SUM_NUJ_0=HLMIJ(NU1,LJ)*RLMIJ(MU1,MJ,LJ)
                ENDIF
C
                SUM_NU1_0=SUM_NU1_0+RLMF_0*SUM_NUJ_0
                SUM_NU1_1=SUM_NU1_1+RLMF_0*SUM_NUJ_1
C
              ENDDO
C
            ENDDO
C
            TAU(INDJP,ILF,JATL)=TAU(INDJP,ILF,JATL)+CONJG(CJ)*TL_J*SUM_N
     &U1_1
            TAU(INDJ,ILF,JATL)=TAU(INDJ,ILF,JATL)+CJ*TL_J*SUM_NU1_0
C
            IF(NPATHP.EQ.0) GOTO 45
C
            FM2=FMAX(JORDP)
            XINT1=CABS(CJ*TL_J*SUM_NU1_0)
            XINT2=CABS(CONJG(CJ)*TL_J*SUM_NU1_1)
            XMAX=AMAX1(XINT1,XINT2,XMAX)
            FMAX(JORDP)=AMAX1(FM2,XINT1,XINT2)
            IF((FMAX(JORDP)-FM2).GT.XCOMP) NPMA(JORDP)=NPATH(JORDP)
            IF((IREF.EQ.1).AND.(JORDP.EQ.NCUT)) THEN
                FREF=FMAX(JORDP)
            ENDIF
 45         CONTINUE
          ENDDO
        ENDDO
C
C  Case MF > 0
C
        CF=ONEC
        DO MF=1,LF
          INDF=ILF+MF
          INDFP=ILF-MF
          CF=CF*EXP_F
C
          DO LJ=0,LMJ
            ILJ=LJ*LJ+LJ+1
            NUJMAX=MIN(LJ,NUMAXJ)
            IF(JORDP.EQ.1) THEN
              NU1MAX=MIN(NU1MAX1,LJ)
            ELSE
              NU1MAX=NU1MAX1
            ENDIF
C
            IF(ISPEED.EQ.1) THEN
              TL_J=COEF*TL(LF,1,1,JE)*TL(LJ,1,JTYP,JE)
            ELSE
              TL_J=COEF*TLT(LF,1,1,JE)*TLT(LJ,1,JTYP,JE)
            ENDIF
C
C    Case MJ = 0
C
            SUM_NU1_0=ZEROC
            SUM_NU1_1=ZEROC
C
            DO NU1=0,NU1MAX
              IF(JORDP.GT.1) THEN
                MU1MAX=MIN(LF-NU1,NO1-NU1-NU1)
              ELSE
                MU1MAX=MIN(LF-NU1,NO1-NU1-NU1,LJ)
              ENDIF
C
              DO MU1=-MU1MAX,MU1MAX
                LAMBDA1=LBD(MU1,NU1)
                AMU1=ABS(MU1)
C
                RLMF_1=HLM01(AMU1+NU1,LF)*RLM01(MU1,-MF,LF)
                RLMF_0=HLM01(AMU1+NU1,LF)*RLM01(MU1,MF,LF)
C
                SUM_NUJ_0=ZEROC
C
                IF(JORDP.GT.1) THEN
                  DO NUJ=0,NUJMAX
                    MUJMAX=MIN(LJ,NOJ-NUJ-NUJ)
C
                    SUM_MUJ_0=ZEROC
C
                    DO MUJ=-MUJMAX,MUJMAX
C
                      LAMBDAJ=LBD(MUJ,NUJ)
C
                      SUM_MUJ_0=SUM_MUJ_0+H(LAMBDAJ,LAMBDA1)*RLMIJ(MUJ,0
     &,LJ)
C
                    ENDDO
C
                    SUM_NUJ_0=SUM_NUJ_0+SUM_MUJ_0*HLMIJ(NUJ,LJ)
C
                  ENDDO
                ELSE
                  SUM_NUJ_0=HLMIJ(NU1,LJ)*RLMIJ(MU1,0,LJ)
                ENDIF
C
                SUM_NU1_0=SUM_NU1_0+RLMF_0*SUM_NUJ_0
                SUM_NU1_1=SUM_NU1_1+RLMF_1*SUM_NUJ_0
C
              ENDDO
C
            ENDDO
C
            TAU(ILJ,INDF,JATL)=TAU(ILJ,INDF,JATL)+CF*TL_J*SUM_NU1_0
            TAU(ILJ,INDFP,JATL)=TAU(ILJ,INDFP,JATL)+CONJG(CF)*TL_J*SUM_N
     &U1_1
C
            IF(NPATHP.EQ.0) GOTO 25
C
            FM2=FMAX(JORDP)
            XINT1=CABS(CF*TL_J*SUM_NU1_0)
            XINT2=CABS(CONJG(CF)*TL_J*SUM_NU1_1)
            XMAX=AMAX1(XINT1,XINT2,XMAX)
            FMAX(JORDP)=AMAX1(FM2,XINT1,XINT2)
            IF((FMAX(JORDP)-FM2).GT.XCOMP) NPMA(JORDP)=NPATH(JORDP)
            IF((IREF.EQ.1).AND.(JORDP.EQ.NCUT)) THEN
                FREF=FMAX(JORDP)
            ENDIF
 25         CONTINUE
C
C    Case MJ > 0
C
            CJ=ONEC
            DO MJ=1,LJ
              INDJ=ILJ+MJ
              INDJP=ILJ-MJ
              CJ=CJ*EXP_J
C
              SUM_NU1_0=ZEROC
              SUM_NU1_1=ZEROC
              SUM_NU1_2=ZEROC
              SUM_NU1_3=ZEROC
C
              DO NU1=0,NU1MAX
                IF(JORDP.GT.1) THEN
                  MU1MAX=MIN(LF-NU1,NO1-NU1-NU1)
                ELSE
                  MU1MAX=MIN(LF-NU1,NO1-NU1-NU1,LJ)
                ENDIF
C
                DO MU1=-MU1MAX,MU1MAX
                  LAMBDA1=LBD(MU1,NU1)
                  AMU1=ABS(MU1)
C
                  RLMF_1=HLM01(AMU1+NU1,LF)*RLM01(MU1,-MF,LF)
                  RLMF_0=HLM01(AMU1+NU1,LF)*RLM01(MU1,MF,LF)
C
                  SUM_NUJ_0=ZEROC
                  SUM_NUJ_1=ZEROC
C
                  IF(JORDP.GT.1) THEN
                    DO NUJ=0,NUJMAX
                      MUJMAX=MIN(LJ,NOJ-NUJ-NUJ)
C
                      SUM_MUJ_0=ZEROC
                      SUM_MUJ_1=ZEROC
C
                      DO MUJ=-MUJMAX,MUJMAX
C
                        LAMBDAJ=LBD(MUJ,NUJ)
C
                        SUM_MUJ_1=SUM_MUJ_1+H(LAMBDAJ,LAMBDA1)*RLMIJ(MUJ
     &,-MJ,LJ)
                        SUM_MUJ_0=SUM_MUJ_0+H(LAMBDAJ,LAMBDA1)*RLMIJ(MUJ
     &,MJ,LJ)
C
                      ENDDO
C
                      SUM_NUJ_0=SUM_NUJ_0+SUM_MUJ_0*HLMIJ(NUJ,LJ)
                      SUM_NUJ_1=SUM_NUJ_1+SUM_MUJ_1*HLMIJ(NUJ,LJ)
C
                    ENDDO
                  ELSE
                    SUM_NUJ_1=HLMIJ(NU1,LJ)*RLMIJ(MU1,-MJ,LJ)
                    SUM_NUJ_0=HLMIJ(NU1,LJ)*RLMIJ(MU1,MJ,LJ)
                  ENDIF
C
                  SUM_NU1_0=SUM_NU1_0+RLMF_0*SUM_NUJ_0
                  SUM_NU1_1=SUM_NU1_1+RLMF_0*SUM_NUJ_1
                  SUM_NU1_2=SUM_NU1_2+RLMF_1*SUM_NUJ_0
                  SUM_NU1_3=SUM_NU1_3+RLMF_1*SUM_NUJ_1
C
                ENDDO
C
              ENDDO
C
              TAU(INDJP,INDF,JATL)=TAU(INDJP,INDF,JATL)+CF*CONJG(CJ)*TL_
     &J*SUM_NU1_1
              TAU(INDJP,INDFP,JATL)=TAU(INDJP,INDFP,JATL)+CONJG(CF*CJ)*T
     &L_J*SUM_NU1_3
              TAU(INDJ,INDF,JATL)=TAU(INDJ,INDF,JATL)+CF*CJ*TL_J*SUM_NU1
     &_0
              TAU(INDJ,INDFP,JATL)=TAU(INDJ,INDFP,JATL)+CONJG(CF)*CJ*TL_
     &J*SUM_NU1_2
C
              IF(NPATHP.EQ.0) GOTO 15
C
              FM2=FMAX(JORDP)
              XINT1=CABS(CF*CJ*TL_J*SUM_NU1_0)
              XINT2=CABS(CF*CONJG(CJ)*TL_J*SUM_NU1_1)
              XINT3=CABS(CONJG(CF)*CJ*TL_J*SUM_NU1_2)
              XINT4=CABS(CONJG(CF*CJ)*TL_J*SUM_NU1_3)
              XMAX=AMAX1(XINT1,XINT2,XINT3,XINT4,XMAX)
              FMAX(JORDP)=AMAX1(FM2,XINT1,XINT2,XINT3,XINT4)
              IF((FMAX(JORDP)-FM2).GT.XCOMP) NPMA(JORDP)=NPATH(JORDP)
              IF((IREF.EQ.1).AND.(JORDP.EQ.NCUT)) THEN
                  FREF=FMAX(JORDP)
              ENDIF
 15           CONTINUE
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
      IF(NPATHP.EQ.0) GOTO 16
      FMIN(JORDP)=AMIN1(FM1,XMAX)
      IF(XMAX.GT.FMN(NPATHP)) THEN
        CALL LOCATE(FMN,NPATHP,XMAX,JMX)
        DO KF=NPATHP,JMX+2,-1
          FMN(KF)=FMN(KF-1)
          JON(KF)=JON(KF-1)
          PATH(KF)=PATH(KF-1)
          DMN(KF)=DMN(KF-1)
          DO KD=1,10
            JPON(KF,KD)=JPON(KF-1,KD)
          ENDDO
        ENDDO
        FMN(JMX+1)=XMAX
        JON(JMX+1)=JORDP
        PATH(JMX+1)=NPATH(JORDP)
        DMN(JMX+1)=D
        DO KD=1,JORDP
          JPON(JMX+1,KD)=JPOS(KD,3)
        ENDDO
      ENDIF
      IF((FMIN(JORDP)-FM1).LT.-XCOMP) NPMI(JORDP)=NPATH(JORDP)
      IF((IPRINT.EQ.3).AND.(IJ.EQ.1)) THEN
        WRITE(IUSCR,1) JORDP,NPATH(JORDP),XMAX,D,(JPOS(KD,3),KD=1,JORDP)
     &
      ENDIF
C
  16  RETURN
C
   1  FORMAT(9X,I2,2X,E12.6,7X,E12.6,1X,F6.3,1X,10(I3,2X))
C
      END
