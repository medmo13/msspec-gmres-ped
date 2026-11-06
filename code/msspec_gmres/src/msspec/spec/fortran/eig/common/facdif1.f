C
C=======================================================================
C
      SUBROUTINE FACDIF1(VKE,RJ,RJK,THRJ,PHIRJ,BETA,GAMMA,L,M,FSPH,JAT,
     &JE,*)
C
C  This routine computes a spherical wave scattering factor
C
C                                        Last modified : 03/04/2006
C
C      INCLUDE 'spec.inc'
C
      USE DIM_MOD
      USE APPROX_MOD
      USE EXPFAC_MOD
      USE TRANS_MOD
      USE TYPCAL_MOD, I2 => IPHI, I3 => IE, I4 => ITHETA, I5 => IMOD, 
     &I6 => IPOL, I7 => I_CP, I8 => I_EXT, I9 => I_TEST

      DIMENSION PLMM(0:100,0:100)
      DIMENSION D(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1)
C
      COMPLEX HLM(0:NO_ST_M,0:NL_M-1),HLN(0:NO_ST_M,0:NL_M-1),FSPH,RHOJ
      COMPLEX HLM1,HLM2,HLM3,HLM4,ALMU,BLMU,SLP,SNU,SMU,VKE
      COMPLEX RHOJK
C
C
      DATA PI/3.141593/
C
      A=1.
      INTER=0
      IF(ITL.EQ.1) VKE=VK(JE)
      RHOJ=VKE*RJ
      RHOJK=VKE*RJK
      HLM1=(1.,0.)
      HLM2=(1.,0.)
      HLM3=(1.,0.)
      HLM4=(1.,0.)
      IEM=1
      CSTH=COS(BETA)
      IF((IFTHET.EQ.0).OR.(THRJ.LT.0.0001)) THEN
        INTER=1
        BLMU=SQRT(4.*PI/FLOAT(2*L+1))*CEXP((0.,-1.)*M*(PHIRJ-PI))
      ENDIF
      CALL PLM(CSTH,PLMM,LMAX(JAT,JE))
      IF(ISPHER.EQ.0) NO1=0
      IF(ISPHER.EQ.1) THEN
        IF(NO.EQ.8) THEN
          NO1=LMAX(JAT,JE)+1
        ELSE
          NO1=NO
        ENDIF
        CALL POLHAN(ISPHER,NO1,LMAX(JAT,JE),RHOJ,HLM)
         IF(IEM.EQ.0) THEN
          HLM4=HLM(0,L)
        ENDIF
        IF(RJK.GT.0.0001) THEN
          NDUM=0
          CALL POLHAN(ISPHER,NDUM,LMAX(JAT,JE),RHOJK,HLN)
        ENDIF
        CALL DJMN(THRJ,D,L)
        A1=ABS(D(0,M,L))
        IF(((A1.LT.0.0001).AND.(IFTHET.EQ.1)).AND.(INTER.EQ.0)) 
     &  RETURN 1
      ENDIF
      MUMAX=MIN0(L,NO1)
      SMU=(0.,0.)
      DO 10 MU=0,MUMAX
        IF(MOD(MU,2).EQ.0) THEN
          B=1.
        ELSE
          B=-1.
          IF(SIN(BETA).LT.0.) THEN
            A=-1.
          ENDIF
        ENDIF
        IF(ISPHER.LE.1) THEN
          ALMU=(1.,0.)
          C=1.
        ENDIF
        IF(ISPHER.EQ.0) GOTO 40
        IF(INTER.EQ.0) BLMU=CMPLX(D(M,0,L))
        IF(MU.GT.0) THEN
          C=B*FLOAT(L+L+1)/EXPF(MU,L)
          ALMU=(D(M,MU,L)*CEXP((0.,-1.)*MU*GAMMA)+B*
     *          CEXP((0.,1.)*MU*GAMMA)*D(M,-MU,L))/BLMU
        ELSE
          C=1.
          ALMU=CMPLX(D(M,0,L))/BLMU
        ENDIF
  40    SNU=(0.,0.)
        NU1=INT(0.5*(NO1-MU)+0.0001)
        NUMAX=MIN0(NU1,L-MU)
        DO 20 NU=0,NUMAX
          SLP=(0.,0.)
          LPMIN=MAX0(MU,NU)
          DO 30 LP=LPMIN,LMAX(JAT,JE)
            IF(ISPHER.EQ.1) THEN
              HLM1=HLM(NU,LP)
              IF(RJK.GT.0.0001) HLM3=HLN(0,LP)
            ENDIF
            SLP=SLP+FLOAT(2*LP+1)*TL(LP,1,JAT,JE)*HLM1*PLMM(LP,
     &      MU)*HLM3
  30      CONTINUE
          IF(ISPHER.EQ.1) THEN
            HLM2=HLM(MU+NU,L)
          ENDIF
          SNU=SNU+SLP*HLM2
  20    CONTINUE
        SMU=SMU+SNU*C*ALMU*A*B
  10  CONTINUE
      FSPH=SMU/(VKE*HLM4)
C
      RETURN
C
      END
