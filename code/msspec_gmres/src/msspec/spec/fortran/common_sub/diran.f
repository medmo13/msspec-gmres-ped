C
C=======================================================================
C
      SUBROUTINE DIRAN(VINT,ECIN,J_EL)
C
C  This subroutine calculates the direction(s) of the analyzer with
C             or without an angular averaging.
C
C             DIRANA is the internal direction
C             ANADIR is the external direction
C
C             J_EL is the type of electron : 1 ---> photoelectron
C                                            2 ---> Auger electron
C
C                                            Last modified : 23/03/2006
C
C
      USE DIRECT_MOD
      USE DIRECT_A_MOD
      USE MOYEN_MOD
      USE MOYEN_A_MOD
      USE OUTUNITS_MOD
      USE TESTS_MOD , I1 => ITEST, I2 => ISORT1, N1 => NPATHP, I3 => ISO
     &M
      USE TYPCAL_MOD
      USE TYPCAL_A_MOD
C
      COMPLEX COEF,IC
      DATA PI,PIS2,PIS180 /3.141593,1.570796,0.017453/
C
      IC=(0.,1.)
C
      IF(J_EL.EQ.1) THEN
        ANADIR(1,1)=SIN(RTHEXT)*COS(RPHI)
        ANADIR(2,1)=SIN(RTHEXT)*SIN(RPHI)
        ANADIR(3,1)=COS(RTHEXT)
        IF((ABS(I_EXT).LE.1).AND.(I_TEST.NE.2)) THEN
          CALL REFRAC(VINT,ECIN,RTHEXT,RTHINT)
        ELSE
          RTHINT=RTHEXT
        ENDIF
        IF((IPRINT.GT.0).AND.(I_EXT.NE.2)) THEN
          DTHEXT=RTHEXT/PIS180
          DTHINT=RTHINT/PIS180
          IF(I_TEST.NE.2) WRITE(IUO1,20) DTHEXT,DTHINT
        ENDIF
        DIRANA(1,1)=SIN(RTHINT)*COS(RPHI)
        DIRANA(2,1)=SIN(RTHINT)*SIN(RPHI)
        DIRANA(3,1)=COS(RTHINT)
        THETAR(1)=RTHINT
        PHIR(1)=RPHI
C
C The change in the definition below is necessary as RPHI is
C   used to define the rotation axis of the direction of the detector
C   when doing polar variations
C
        IF(ITHETA.EQ.1) THEN
          IF(RPHI.GT.PIS2) THEN
            RTHEXT=-RTHEXT
            RPHI=RPHI-PI
          ELSEIF(RPHI.LT.-PIS2) THEN
            RTHEXT=-RTHEXT
            RPHI=RPHI+PI
          ENDIF
        ENDIF
C
        IF(IMOY.GE.1) THEN
          N=2**(IMOY-1)
          S=SIN(ACCEPT*PI/180.)
          RN=FLOAT(N)
          J=1
          DO K1=-N,N
            RK1=FLOAT(K1)
            DO K2=-N,N
              RK2=FLOAT(K2)
              D=SQRT(RK1*RK1+RK2*RK2)
              IF((D-RN).GT.0.000001) GOTO 10
              IF((K1.EQ.0).AND.(K2.EQ.0)) GOTO 10
              C=SQRT(RN*RN-(RK1*RK1+RK2*RK2)*S*S)
              J=J+1
C
              ANADIR(1,J)=(RK1*S*COS(RTHEXT)*COS(RPHI) -RK2*S*SIN(RPHI)+
     &C*ANADIR(1,1))/RN
              ANADIR(2,J)=(RK1*S*COS(RTHEXT)*SIN(RPHI) +RK2*S*COS(RPHI)+
     &C*ANADIR(2,1))/RN
              ANADIR(3,J)=(-RK1*S*SIN(RTHEXT) +C*ANADIR(3,1))/RN
              THETA_R=ACOS(ANADIR(3,J))
              COEF=ANADIR(1,J)+IC*ANADIR(2,J)
              CALL ARCSIN(COEF,ANADIR(3,J),PHI_R)
              IF((ABS(I_EXT).LE.1).AND.(I_TEST.NE.2)) THEN
                CALL REFRAC(VINT,ECIN,THETA_R,THINT_R)
              ELSE
                THINT_R=THETA_R
              ENDIF
C
              DIRANA(1,J)=SIN(THINT_R)*COS(PHI_R)
              DIRANA(2,J)=SIN(THINT_R)*SIN(PHI_R)
              DIRANA(3,J)=COS(THINT_R)
C
              THETAR(J)=THINT_R
              PHIR(J)=PHI_R
  10          CONTINUE
            ENDDO
          ENDDO
        ENDIF
C
      ELSEIF(J_EL.EQ.2) THEN
        ANADIR_A(1,1)=SIN(RTHEXT_A)*COS(RPHI_A)
        ANADIR_A(2,1)=SIN(RTHEXT_A)*SIN(RPHI_A)
        ANADIR_A(3,1)=COS(RTHEXT_A)
        IF((ABS(I_EXT_A).LE.1).AND.(I_TEST_A.NE.2)) THEN
          CALL REFRAC(VINT,ECIN,RTHEXT_A,RTHINT_A)
        ELSE
          RTHINT_A=RTHEXT_A
        ENDIF
        IF((IPRINT.GT.0).AND.(I_EXT_A.NE.2)) THEN
          DTHEXT_A=RTHEXT_A/PIS180
          DTHINT_A=RTHINT_A/PIS180
          IF(I_TEST_A.NE.2) WRITE(IUO1,21) DTHEXT_A,DTHINT_A
        ENDIF
        DIRANA_A(1,1)=SIN(RTHINT_A)*COS(RPHI_A)
        DIRANA_A(2,1)=SIN(RTHINT_A)*SIN(RPHI_A)
        DIRANA_A(3,1)=COS(RTHINT_A)
        THETAR_A(1)=RTHINT_A
        PHIR_A(1)=RPHI_A
C
C The change in the definition below is necessary as RPHI is
C   used to define the rotation axis of the direction of the detector
C   when doing polar variations
C
        IF(ITHETA_A.EQ.1) THEN
          IF(RPHI_A.GT.PIS2) THEN
            RTHEXT_A=-RTHEXT_A
            RPHI_A=RPHI_A-PI
          ELSEIF(RPHI_A.LT.-PIS2) THEN
            RTHEXT_A=-RTHEXT_A
            RPHI_A=RPHI_A+PI
          ENDIF
        ENDIF
C
        IF(IMOY_A.GE.1) THEN
          N=2**(IMOY_A-1)
          S=SIN(ACCEPT_A*PI/180.)
          RN=FLOAT(N)
          J=1
          DO K1=-N,N
            RK1=FLOAT(K1)
            DO K2=-N,N
              RK2=FLOAT(K2)
              D=SQRT(RK1*RK1+RK2*RK2)
              IF((D-RN).GT.0.000001) GOTO 15
              IF((K1.EQ.0).AND.(K2.EQ.0)) GOTO 15
              C=SQRT(RN*RN-(RK1*RK1+RK2*RK2)*S*S)
              J=J+1
C
              ANADIR_A(1,J)=(RK1*S*COS(RTHEXT_A)*COS(RPHI_A) -RK2*S*SIN(
     &RPHI_A)+C*ANADIR_A(1,1))/RN
              ANADIR_A(2,J)=(RK1*S*COS(RTHEXT_A)*SIN(RPHI_A) +RK2*S*COS(
     &RPHI_A)+C*ANADIR_A(2,1))/RN
              ANADIR_A(3,J)=(-RK1*S*SIN(RTHEXT_A) +C*ANADIR_A(3,1))/RN
              THETA_R_A=ACOS(ANADIR_A(3,J))
              COEF=ANADIR_A(1,J)+IC*ANADIR_A(2,J)
              CALL ARCSIN(COEF,ANADIR_A(3,J),PHI_R_A)
              IF((ABS(I_EXT_A).LE.1).AND.(I_TEST_A.NE.2)) THEN
                CALL REFRAC(VINT,ECIN,THETA_R_A,THINT_R_A)
              ELSE
                THINT_R_A=THETA_R_A
              ENDIF
C
              DIRANA_A(1,J)=SIN(THINT_R_A)*COS(PHI_R_A)
              DIRANA_A(2,J)=SIN(THINT_R_A)*SIN(PHI_R_A)
              DIRANA_A(3,J)=COS(THINT_R_A)
C
              THETAR_A(J)=THINT_R_A
              PHIR_A(J)=PHI_R_A
  15          CONTINUE
            ENDDO
          ENDDO
        ENDIF
C
      ENDIF
C
  20  FORMAT(/,10X,'PHOTOELECTRON EXTERNAL THETA  =',F7.2,5X,'INTERNAL T
     &HETA =', F7.2)
  21  FORMAT(/,10X,'AUGER ELECTRON EXTERNAL THETA =',F7.2,5X,'INTERNAL T
     &HETA =', F7.2)
C
      RETURN
C
      END
