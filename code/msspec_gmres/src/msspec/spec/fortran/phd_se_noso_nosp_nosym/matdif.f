C
C=======================================================================
C
      SUBROUTINE MATDIF(NO,ND,LF,JTYP,KTYP,JE,I_ABS,ISPEED,ISPHER,A21,B2
     &1,C21,RHO1,RHO2)
C
C  This routine calculates the Rehr-Albers scattering matrix
C     F_{LAMBDA1,LAMBDA2}. The result is stored in the COMMON block
C     /SCATMAT/ as F21(NSPIN2_M,NLAMBDA_M,NLAMBDA_M,NDIF_M).
C
C                                         Last modified :  3 Aug 2007
C
      USE DIM_MOD
C
      USE EXPFAC_MOD
      USE LBD_MOD
      USE LINLBD_MOD
      USE RA_MOD
      USE SCATMAT_MOD
      USE TRANS_MOD
      USE TLDW_MOD
C
      REAL RLM(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1)
C
      COMPLEX HLM1(0:NO_ST_M,0:NL_M-1),HLM2(0:NO_ST_M,0:NL_M-1)
      COMPLEX SL,RHO1,RHO2,IC,ZEROC,ONEC,ONEOVK
      COMPLEX SL_2_1,SL_2_2
      COMPLEX EXP1,EXP2,PROD1,PROD2
C
      DATA PI,SMALL /3.141593,0.0001/
C
      IC=(0.,1.)
      ZEROC=(0.,0.)
      ONEC=(1.,0.)
      ONEOVK=1./VK(JE)
      IB=0
      LMJ=LMAX(JTYP,JE)
      IF(ABS(ABS(B21)-PI).LT.SMALL) IB=-1
      IF(ABS(B21).LT.SMALL) IB=1
      IF(NO.EQ.8) THEN
        NN2=LMAX(JTYP,JE)+1
      ELSE
        NN2=NO
      ENDIF
C
C  NO is atom-dependent and is decreased with the rank of the scatterer
C    in the path when I_NO > 0. Here LAMBDA1 depends on the scatterer JTYP
C    while LAMBDA2 depends on the next atom (KTYP) in the path
C
      IF(I_NO.EQ.0) THEN
        NO1=N_RA(JTYP)
        NO2=N_RA(KTYP)
      ELSE
        NO1=MAX(N_RA(JTYP)-(ND-1)/I_NO,0)
        NO2=MAX(N_RA(KTYP)-ND/I_NO,0)
      ENDIF
      IF(I_ABS.EQ.0) THEN
        NUMAX1=NO1/2
        NUMAX2=NO2/2
      ELSEIF(I_ABS.EQ.1) THEN
        NUMAX1=MIN0(LF,NO1/2)
        NUMAX2=NO2/2
      ELSEIF(I_ABS.EQ.2) THEN
        NUMAX1=NO1/2
        NUMAX2=MIN0(LF,NO2/2)
      ENDIF
      LBDM(1,ND)=(NO1+1)*(NO1+2)/2
      LBDM(2,ND)=(NO2+1)*(NO2+2)/2
C
      EXP2=-EXP(-IC*A21)
      EXP1=EXP(-IC*C21)
C
      DO LAMBDA1=1,LBDMAX
        DO LAMBDA2=1,LBDMAX
          F21(1,LAMBDA2,LAMBDA1,ND)=ZEROC
        ENDDO
      ENDDO
C
      IF(ABS(RHO1-RHO2).GT.SMALL) THEN
         CALL POLHAN(ISPHER,NUMAX1,LMJ,RHO1,HLM1)
         CALL POLHAN(ISPHER,NN2,LMJ,RHO2,HLM2)
         NEQUAL=0
      ELSE
         CALL POLHAN(ISPHER,NN2,LMJ,RHO1,HLM1)
         NEQUAL=1
      ENDIF
C
C  Calculation of the scattering matrix when the scattering angle
C              is different from 0 and pi
C
      IF(IB.EQ.0) THEN
        CALL DJMN(B21,RLM,LMJ)
        DO NU1=0,NUMAX1
          MUMAX1=NO1-2*NU1
          IF(I_ABS.EQ.1) MUMAX1=MIN(LF-NU1,MUMAX1)
          DO NU2=0,NUMAX2
            MUMAX2=NO2-2*NU2
C
C  Case MU1 = 0
C
            LAMBDA1=LBD(0,NU1)
C
C    Case MU2 = 0
C
            LAMBDA2=LBD(0,NU2)
            LMIN=MAX(NU1,NU2)
            SL=ZEROC
            DO L=LMIN,LMJ
              IF(NEQUAL.EQ.1) THEN
                HLM2(NU2,L)=HLM1(NU2,L)
              ENDIF
              IF(ISPEED.EQ.1) THEN
                SL=SL+FLOAT(L+L+1)*RLM(0,0,L)*TL(L,1,JTYP,JE)*HLM1(NU1,L
     &)*HLM2(NU2,L)
              ELSE
                SL=SL+FLOAT(L+L+1)*RLM(0,0,L)*TLT(L,1,JTYP,JE)*HLM1(NU1,
     &L)*HLM2(NU2,L)
              ENDIF
            ENDDO
            F21(1,LAMBDA2,LAMBDA1,ND)=SL*ONEOVK
C
C    Case MU2 > 0
C
            PROD2=ONEC
            SIG2=1.
            DO MU2=1,MUMAX2
              LAMBDA2_1=LBD(MU2,NU2)
              LAMBDA2_2=LBD(-MU2,NU2)
              PROD2=PROD2*EXP2
              SIG2=-SIG2
              LMIN=MAX(NU1,MU2+NU2)
              SL=ZEROC
              DO L=LMIN,LMJ
                IF(NEQUAL.EQ.1) THEN
                  HLM2(MU2+NU2,L)=HLM1(MU2+NU2,L)
                ENDIF
                C1=EXPF(0,L)/EXPF(MU2,L)
                IF(ISPEED.EQ.1) THEN
                  SL=SL+FLOAT(L+L+1)*RLM(MU2,0,L)*C1*TL(L,1,JTYP,JE)*HLM
     &1(NU1,L)*HLM2(MU2+NU2,L)
                ELSE
                  SL=SL+FLOAT(L+L+1)*RLM(MU2,0,L)*C1*TLT(L,1,JTYP,JE)*HL
     &M1(NU1,L)*HLM2(MU2+NU2,L)
                ENDIF
              ENDDO
              F21(1,LAMBDA2_1,LAMBDA1,ND)=SL*PROD2*ONEOVK*SIG2
              F21(1,LAMBDA2_2,LAMBDA1,ND)=SL*ONEOVK/PROD2
            ENDDO
C
C  Case MU1 > 0
C
            PROD1=ONEC
            SIG1=1.
            DO MU1=1,MUMAX1
              LAMBDA1_1=LBD(MU1,NU1)
              LAMBDA1_2=LBD(-MU1,NU1)
              PROD1=PROD1*EXP1
              SIG1=-SIG1
C
C    Case MU2 = 0
C
              LAMBDA2=LBD(0,NU2)
              LMIN=MAX(MU1,NU1,NU2)
              SL=ZEROC
              DO L=LMIN,LMJ
                IF(NEQUAL.EQ.1) THEN
                  HLM2(NU2,L)=HLM1(NU2,L)
                ENDIF
                C1=EXPF(MU1,L)/EXPF(0,L)
                IF(ISPEED.EQ.1) THEN
                  SL=SL+FLOAT(L+L+1)*RLM(0,MU1,L)*C1*TL(L,1,JTYP,JE)*HLM
     &1(NU1,L)*HLM2(NU2,L)
                ELSE
                  SL=SL+FLOAT(L+L+1)*RLM(0,MU1,L)*C1*TLT(L,1,JTYP,JE)*HL
     &M1(NU1,L)*HLM2(NU2,L)
                ENDIF
              ENDDO
              F21(1,LAMBDA2,LAMBDA1_1,ND)=SL*PROD1*ONEOVK*SIG1
              F21(1,LAMBDA2,LAMBDA1_2,ND)=SL*ONEOVK/PROD1
C
C    Case MU2 > 0
C
              PROD2=ONEC
              SIG2=SIG1
              DO MU2=1,MUMAX2
                LAMBDA2_1=LBD(MU2,NU2)
                LAMBDA2_2=LBD(-MU2,NU2)
                PROD2=PROD2*EXP2
                SIG2=-SIG2
                LMIN=MAX(MU1,NU1,MU2+NU2)
                SL_2_1=ZEROC
                SL_2_2=ZEROC
                DO L=LMIN,LMJ
                  IF(NEQUAL.EQ.1) THEN
                    HLM2(MU2+NU2,L)=HLM1(MU2+NU2,L)
                  ENDIF
                  C1=EXPF(MU1,L)/EXPF(MU2,L)
                  IF(ISPEED.EQ.1) THEN
                    SL=FLOAT(L+L+1)*C1*TL(L,1,JTYP,JE)*HLM1(NU1,L)*HLM2(
     &MU2+NU2,L)
                  ELSE
                    SL=FLOAT(L+L+1)*C1*TLT(L,1,JTYP,JE)*HLM1(NU1,L)*HLM2
     &(MU2+NU2,L)
                  ENDIF
                  SL_2_1=SL_2_1+SL*RLM(MU2,-MU1,L)
                  SL_2_2=SL_2_2+SL*RLM(MU2,MU1,L)
                ENDDO
                F21(1,LAMBDA2_1,LAMBDA1_1,ND)=SL_2_2*PROD1*PROD2*ONEOVK*
     &SIG2
                F21(1,LAMBDA2_2,LAMBDA1_1,ND)=SL_2_1*PROD1*ONEOVK/PROD2
                F21(1,LAMBDA2_1,LAMBDA1_2,ND)=SL_2_1*ONEOVK*PROD2*SIG2/P
     &ROD1
                F21(1,LAMBDA2_2,LAMBDA1_2,ND)=SL_2_2*ONEOVK/(PROD1*PROD2
     &)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C
C  Calculation of the scattering matrix when the scattering angle
C     is equal to 0 (forward scattering) or pi (backscattering)
C
      ELSEIF(IB.EQ.1) THEN
        DO NU1=0,NUMAX1
          DO NU2=0,NUMAX2
            MUMAX1=MIN0(NO1-2*NU1,NO1-2*NU2)
            IF(I_ABS.EQ.1) MUMAX1=MIN0(LF-NU1,MUMAX1)
C
C  Case MU = 0
C
            LAMBDA1=LBD(0,NU1)
            LAMBDA2=LBD(0,NU2)
            LMIN=MAX(NU1,NU2)
            SL=ZEROC
            DO L=LMIN,LMJ
              IF(NEQUAL.EQ.1) THEN
                 HLM2(NU2,L)=HLM1(NU2,L)
              ENDIF
              IF(ISPEED.EQ.1) THEN
                SL=SL+FLOAT(L+L+1)*TL(L,1,JTYP,JE)*HLM1(NU1,L)*HLM2(NU2,
     &L)
              ELSE
                SL=SL+FLOAT(L+L+1)*TLT(L,1,JTYP,JE)*HLM1(NU1,L)*HLM2(NU2
     &,L)
              ENDIF
            ENDDO
            F21(1,LAMBDA2,LAMBDA1,ND)=SL*ONEOVK
C
C  Case MU > 0
C
            CST1=1.
            DO MU=1,MUMAX1
              LAMBDA1=LBD(MU,NU2)
              LAMBDA2=LBD(-MU,NU2)
              CST1=-CST1
              LMIN=MAX(NU1,MU+NU2)
              SL=ZEROC
              DO L=LMIN,LMJ
                IF(NEQUAL.EQ.1) THEN
                   HLM2(MU+NU2,L)=HLM1(MU+NU2,L)
                ENDIF
                IF(ISPEED.EQ.1) THEN
                  SL=SL+FLOAT(L+L+1)*CST1*TL(L,1,JTYP,JE)*HLM1(NU1,L)*HL
     &M2(MU+NU2,L)
                ELSE
                  SL=SL+FLOAT(L+L+1)*CST1*TLT(L,1,JTYP,JE)*HLM1(NU1,L)*H
     &LM2(MU+NU2,L)
                ENDIF
              ENDDO
              F21(1,LAMBDA1,LAMBDA1,ND)=SL*ONEOVK
              F21(1,LAMBDA2,LAMBDA2,ND)=SL*ONEOVK
            ENDDO
          ENDDO
        ENDDO
      ELSEIF(IB.EQ.-1) THEN
        DO NU1=0,NUMAX1
          DO NU2=0,NUMAX2
            MUMAX1=MIN0(NO1-2*NU1,NO1-2*NU2)
            IF(I_ABS.EQ.1) MUMAX1=MIN0(LF-NU1,MUMAX1)
C
C  Case MU = 0
C
            LAMBDA1=LBD(0,NU1)
            LAMBDA2=LBD(0,NU2)
            LMIN=MAX(NU1,NU2)
            SL=ZEROC
            DO L=LMIN,LMJ
              IF(NEQUAL.EQ.1) THEN
                 HLM2(NU2,L)=HLM1(NU2,L)
              ENDIF
              IF(MOD(L,2).EQ.0) THEN
                CST2=1.0
              ELSE
                CST2=-1.0
              ENDIF
              IF(ISPEED.EQ.1) THEN
                SL=SL+FLOAT(L+L+1)*CST2*TL(L,1,JTYP,JE)*HLM1(NU1,L)*HLM2
     &(NU2,L)
              ELSE
                SL=SL+FLOAT(L+L+1)*CST2*TLT(L,1,JTYP,JE)*HLM1(NU1,L)*HLM
     &2(NU2,L)
              ENDIF
            ENDDO
            F21(1,LAMBDA2,LAMBDA1,ND)=SL*ONEOVK
C
C  Case MU > 0
C
            CST1=1.
            DO MU=1,MUMAX1
              MUP=-MU
              LAMBDA1_1=LBD(MUP,NU1)
              LAMBDA1_2=LBD(-MUP,NU1)
              LAMBDA2_1=LBD(MU,NU2)
              LAMBDA2_2=LBD(-MU,NU2)
              CST1=-CST1
              LMIN=MAX(NU1,MU+NU2)
              SL=ZEROC
              DO L=LMIN,LMJ
                IF(NEQUAL.EQ.1) THEN
                   HLM2(MU+NU2,L)=HLM1(MU+NU2,L)
                ENDIF
                IF(MOD(L,2).EQ.0) THEN
                  CST2=CST1
                ELSE
                  CST2=-CST1
                ENDIF
                IF(ISPEED.EQ.1) THEN
                  SL=SL+FLOAT(L+L+1)*CST2*TL(L,1,JTYP,JE)*HLM1(NU1,L)*HL
     &M2(MU+NU2,L)
                ELSE
                  SL=SL+FLOAT(L+L+1)*CST2*TLT(L,1,JTYP,JE)*HLM1(NU1,L)*H
     &LM2(MU+NU2,L)
                ENDIF
              ENDDO
              F21(1,LAMBDA2_1,LAMBDA1_1,ND)=SL*ONEOVK
              F21(1,LAMBDA2_2,LAMBDA1_2,ND)=SL*ONEOVK
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
C
      END
