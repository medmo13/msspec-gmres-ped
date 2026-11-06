C
C=======================================================================
C
      SUBROUTINE COUMAT(ITL,MI,LF,MF,DELTA,RADIAL,MATRIX)
C
C  This routine calculates the spin-independent PhD optical matrix
C       elements for dipolar excitations. It is stored in
C                          MATRIX(JDIR,JPOL)
C
C  Here, the conventions are :
C
C             IPOL=1  : linearly polarized light
C             IPOL=2  : circularly polarized light
C
C             JPOL=1  : +/x polarization for circular/linear light
C             JPOL=2  : -/y polarization for circular/linear light
C
C  When IDICHR=0, JDIR = 1,2 and 3 correspond respectively to the x,y
C       and z directions for the linear polarization. But for IDICHR=1,
C       these basis directions are those of the position of the light.
C
C                                     Last modified :  8 Dec 2008
C
      USE DIM_MOD
C
      USE INIT_L_MOD , L2 => NNL, L3 => LF1, L4 => LF2, L5 => ISTEP_LF
      USE SPIN_MOD , I1 => ISPIN, N1 => NSPIN, N2 => NSPIN2, I2 => ISFLI
     &P, I8 => IR_DIA, N3 => NSTEP
      USE TYPCAL_MOD , I3 => IPHI, I4 => IE, I5 => ITHETA, I6 => IFTHET,
     & I7 => IMOD, I9 => I_CP, I10 => I_EXT
C
      COMPLEX MATRIX(3,2),SUM_1,SUM_2,DELTA,YLM(3,-1:1),RADIAL
      COMPLEX ONEC,IC,IL,COEF,PROD
C
      REAL RLM(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1),GNT(0:N_GAUNT)
      REAL THETA(3),PHI(3)
C
      DATA PI4S3,C_LIN,SQR2 /4.188790,1.447202,1.414214/
      DATA PIS2 /1.570796/
C
      ONEC=(1.,0.)
      IC=(0.,1.)
C
      IF(INITL.EQ.0) GOTO 2
C
      M=MF-MI
C
      IF(MOD(LF,4).EQ.0) THEN
        IL=ONEC
      ELSEIF(MOD(LF,4).EQ.1) THEN
        IL=IC
      ELSEIF(MOD(LF,4).EQ.2) THEN
        IL=-ONEC
      ELSEIF(MOD(LF,4).EQ.3) THEN
        IL=-IC
      ENDIF
C
      CALL GAUNT(LI,MI,LF,MF,GNT)
C
      IF(ITL.EQ.0) THEN
c        COEF=CEXP(IC*DELTA)*CONJG(IL)
        COEF=CEXP(IC*DELTA)*IL
      ELSE
        IF(IDICHR.EQ.0) THEN
c          COEF=PI4S3*CONJG(IL)
          COEF=PI4S3*IL
        ELSE
c          COEF=C_LIN*CONJG(IL)
          COEF=C_LIN*IL
        ENDIF
      ENDIF
C
      PROD=COEF*RADIAL*GNT(1)
C
      IF(IDICHR.EQ.0) THEN
        YLM(1,-1)=(0.345494,0.)
        YLM(1,0)=(0.,0.)
        YLM(1,1)=(-0.345494,0.)
        YLM(2,-1)=(0.,-0.345494)
        YLM(2,0)=(0.,0.)
        YLM(2,1)=(0.,-0.345494)
        YLM(3,-1)=(0.,0.)
        YLM(3,0)=(0.488602,0.)
        YLM(3,1)=(0.,0.)
C
        DO JDIR=1,3
          MATRIX(JDIR,1)=PROD*CONJG(YLM(JDIR,M))
        ENDDO
C
      ELSEIF(IDICHR.GE.1) THEN
C
        THETA(1)=PIS2
        PHI(1)=0.
        THETA(2)=PIS2
        PHI(2)=PIS2
        THETA(3)=0.
        PHI(3)=0.
C
        DO JDIR=1,3
          CALL DJMN(THETA(JDIR),RLM,1)
          SUM_1=RLM(-1,M,1)*PROD*CEXP((0.,-1.)*M*PHI(JDIR))
          SUM_2=RLM(1,M,1)*PROD*CEXP((0.,-1.)*M*PHI(JDIR))
          IF(IPOL.EQ.2) THEN
            MATRIX(JDIR,1)=SQR2*SUM_1
            MATRIX(JDIR,2)=SQR2*SUM_2
          ELSEIF(ABS(IPOL).EQ.1) THEN
            MATRIX(JDIR,1)=(SUM_2-SUM_1)
            MATRIX(JDIR,2)=(SUM_2+SUM_1)*IC
          ENDIF
        ENDDO
      ENDIF
      GOTO 1
C
   2  DO JDIR=1,3
        MATRIX(JDIR,1)=ONEC
        MATRIX(JDIR,2)=ONEC
      ENDDO
C
   1  RETURN
C
      END
