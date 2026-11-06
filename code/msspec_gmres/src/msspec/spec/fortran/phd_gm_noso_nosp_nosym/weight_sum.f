C
C=======================================================================
C
      SUBROUTINE WEIGHT_SUM(ISOM,I_EXT,I_EXT_A,JEL)
C
C   This subroutine performs a weighted sum of the results
C     corresponding to different directions of the detector.
C     The directions and weights are read from an external input file
C
C   JEL is the electron undetected (i.e. for which the outgoing
C     directions are integrated over the unit sphere). It is always
C     1 for one electron spectroscopies (PHD). For APECS, It can be
C     1 (photoelectron) or 2 (Auger electron) or even 0 (no electron
C     detected)
C
C                                         Last modified : 31 Jan 2007
C
      USE DIM_MOD
      USE INFILES_MOD
      USE INUNITS_MOD
      USE OUTUNITS_MOD
C
C
      PARAMETER(N_MAX=5810,NPM=20)
C
      REAL*4 W(N_MAX),W_A(N_MAX),ECIN(NE_M)
      REAL*4 DTHETA(N_MAX),DPHI(N_MAX),DTHETAA(N_MAX),DPHIA(N_MAX)
      REAL*4 SR_1,SF_1,SR_2,SF_2
      REAL*4 SUMR_1(NPM,NE_M,N_MAX),SUMR_2(NPM,NE_M,N_MAX)
      REAL*4 SUMF_1(NPM,NE_M,N_MAX),SUMF_2(NPM,NE_M,N_MAX)
C
      CHARACTER*3 SPECTRO,SPECTRO2
      CHARACTER*5 LIKE
      CHARACTER*13 OUTDATA
C
C
C
C
      DATA JVOL,JTOT/0,-1/
      DATA LIKE /'-like'/
C
      REWIND IUO2
C
      READ(IUO2,15) SPECTRO,OUTDATA
      IF(SPECTRO.NE.'APC') THEN
        READ(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        READ(IUO2,8) NPHI,NTHETA,NE,NPLAN,ISOM
        SPECTRO2='XAS'
      ELSE
        READ(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        READ(IUO2,9) ISPIN_A,IDICHR_A,I_SO_A,ISFLIP_A,ICHKDIR_A,IPHI_A,I
     &THETA_A,IE_A
        READ(IUO2,8) NPHI,NTHETA,NE,NPLAN,ISOM
        READ(IUO2,8) NPHI_A,NTHETA_A
        IF(JEL.EQ.1) THEN
          SPECTRO2='AED'
        ELSEIF(JEL.EQ.2) THEN
          SPECTRO2='PHD'
        ELSEIF(JEL.EQ.0) THEN
          SPECTRO2='XAS'
        ENDIF
      ENDIF
C
      IF(NPLAN.GT.NPM) THEN
        WRITE(IUO1,4) NPLAN+2
        STOP
      ENDIF
C
C  Reading the number of angular points
C
      IF(SPECTRO.NE.'APC') THEN
        OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
        READ(IUI6,1) N_POINTS
        READ(IUI6,5) I_DIM,N_DUM1,N_DUM2
        N_POINTS_A=1
      ELSE
        IF(JEL.EQ.1) THEN
          OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
          READ(IUI6,1) N_POINTS
          READ(IUI6,5) I_DIM,N_DUM1,N_DUM2
          IF(I_EXT_A.EQ.0) THEN
            N_POINTS_A=NTHETA_A*NPHI_A
          ELSE
            OPEN(UNIT=IUI9, FILE=INFILE9, STATUS='OLD')
            READ(IUI9,1) N_POINTS_A
            READ(IUI9,5) I_DIM,N_DUM1,N_DUM2
          ENDIF
          NTHETA0=NTHETA_A
          NPHI0=NPHI_A
        ELSEIF(JEL.EQ.2) THEN
          OPEN(UNIT=IUI9, FILE=INFILE9, STATUS='OLD')
          READ(IUI9,1) N_POINTS_A
          READ(IUI9,5) I_DIM,N_DUM1,N_DUM2
          IF(I_EXT.EQ.0) THEN
            N_POINTS=NTHETA*NPHI
          ELSE
            OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
            READ(IUI6,1) N_POINTS
            READ(IUI6,5) I_DIM,N_DUM1,N_DUM2
          ENDIF
          NTHETA0=NTHETA
          NPHI0=NPHI
        ELSEIF(JEL.EQ.0) THEN
          OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
          OPEN(UNIT=IUI9, FILE=INFILE9, STATUS='OLD')
          READ(IUI6,1) N_POINTS
          READ(IUI9,1) N_POINTS_A
          READ(IUI6,5) I_DIM,N_DUM1,N_DUM2
          READ(IUI9,5) I_DIM,N_DUM1,N_DUM2
        ENDIF
      ENDIF
C
      IF(SPECTRO.NE.'APC') THEN
        NANGLE=1
      ELSE
        IF(JEL.EQ.1) THEN
          NANGLE=N_POINTS_A
        ELSEIF(JEL.EQ.2) THEN
          NANGLE=N_POINTS
        ELSEIF(JEL.EQ.0) THEN
          NANGLE=1
        ENDIF
      ENDIF
C
C  Initialization of the arrays
C
      DO JE=1,NE
        DO JANGLE=1,NANGLE
          DO JPLAN=1,NPLAN+2
            SUMR_1(JPLAN,JE,JANGLE)=0.
            SUMF_1(JPLAN,JE,JANGLE)=0.
            IF(IDICHR.GT.0) THEN
              SUMR_2(JPLAN,JE,JANGLE)=0.
              SUMF_2(JPLAN,JE,JANGLE)=0.
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
C  Reading of the data to be angle integrated
C
      DO JE=1,NE
C
        DO JANGLE=1,N_POINTS
          IF(I_EXT.NE.0) READ(IUI6,2) TH,PH,W(JANGLE)
          DO JANGLE_A=1,N_POINTS_A
            IF((I_EXT_A.NE.0).AND.(JANGLE.EQ.1)) THEN
               READ(IUI9,2) THA,PHA,W_A(JANGLE_A)
            ENDIF
C
            DO JPLAN=1,NPLAN+2
C
              IF(IDICHR.EQ.0) THEN
                IF(SPECTRO.NE.'APC') THEN
                  READ(IUO2,3) JDUM,DTHETA(JANGLE),DPHI(JANGLE),ECIN(JE)
     &,SR_1,SF_1
                ELSE
                  READ(IUO2,13) JDUM,DTHETA(JANGLE),DPHI(JANGLE),ECIN(JE
     &),DTHETAA(JANGLE_A),DPHIA(JANGLE_A),SR_1,SF_1
                ENDIF
              ELSE
                IF(SPECTRO.NE.'APC') THEN
                  READ(IUO2,23) JDUM,DTHETA(JANGLE),DPHI(JANGLE),ECIN(JE
     &),SR_1,SF_1,SR_2,SF_2
                ELSE
                  READ(IUO2,24) JDUM,DTHETA(JANGLE),DPHI(JANGLE),ECIN(JE
     &),DTHETAA(JANGLE_A),DPHIA(JANGLE_A),SR_1,SF_1,SR_2,SF_2
                ENDIF
              ENDIF
C
              IF(JEL.EQ.1) THEN
                SUMR_1(JPLAN,JE,JANGLE_A)=SUMR_1(JPLAN,JE,JANGLE_A)+SR_1
     &*W(JANGLE)
                SUMF_1(JPLAN,JE,JANGLE_A)=SUMF_1(JPLAN,JE,JANGLE_A)+SF_1
     &*W(JANGLE)
              ELSEIF(JEL.EQ.2) THEN
                SUMR_1(JPLAN,JE,JANGLE)=SUMR_1(JPLAN,JE,JANGLE)+SR_1*W_A
     &(JANGLE_A)
                SUMF_1(JPLAN,JE,JANGLE)=SUMF_1(JPLAN,JE,JANGLE)+SF_1*W_A
     &(JANGLE_A)
              ELSEIF(JEL.EQ.0) THEN
                SUMR_1(JPLAN,JE,1)=SUMR_1(JPLAN,JE,1)+SR_1*W(JANGLE)*W_A
     &(JANGLE_A)
                SUMF_1(JPLAN,JE,1)=SUMF_1(JPLAN,JE,1)+SF_1*W(JANGLE)*W_A
     &(JANGLE_A)
              ENDIF
              IF(IDICHR.GT.0) THEN
                IF(JEL.EQ.1) THEN
                  SUMR_2(JPLAN,JE,JANGLE_A)=SUMR_2(JPLAN,JE,JANGLE_A)+SR
     &_2*W(JANGLE)
                  SUMF_2(JPLAN,JE,JANGLE_A)=SUMF_2(JPLAN,JE,JANGLE_A)+SF
     &_2*W(JANGLE)
                ELSEIF(JEL.EQ.2) THEN
                  SUMR_2(JPLAN,JE,JANGLE)=SUMR_2(JPLAN,JE,JANGLE)+SR_2*W
     &_A(JANGLE_A)
                  SUMF_2(JPLAN,JE,JANGLE)=SUMF_2(JPLAN,JE,JANGLE)+SF_2*W
     &_A(JANGLE_A)
                ELSEIF(JEL.EQ.0) THEN
                  SUMR_2(JPLAN,JE,1)=SUMR_2(JPLAN,JE,1)+SR_2*W(JANGLE)*W
     &_A(JANGLE_A)
                  SUMF_2(JPLAN,JE,1)=SUMF_2(JPLAN,JE,1)+SF_2*W(JANGLE)*W
     &_A(JANGLE_A)
                ENDIF
              ENDIF
C
            ENDDO
C
          ENDDO
          IF(I_EXT_A.NE.0) THEN
            REWIND IUI9
            READ(IUI9,1) NDUM
            READ(IUI9,1) NDUM
          ENDIF
        ENDDO
C
        IF(I_EXT.NE.0) THEN
          REWIND IUI6
          READ(IUI6,1) NDUM
          READ(IUI6,1) NDUM
        ENDIF
      ENDDO
C
      CLOSE(IUI6)
      CLOSE(IUI9)
      REWIND IUO2
C
      WRITE(IUO2,16) SPECTRO2,LIKE,SPECTRO,OUTDATA
      IF((SPECTRO.NE.'APC').OR.(JEL.EQ.0)) THEN
        WRITE(IUO2,19) ISPIN,IDICHR,I_SO,ISFLIP
        WRITE(IUO2,18) NE,NPLAN,ISOM
      ELSEIF(JEL.EQ.1) THEN
        WRITE(IUO2,20) ISPIN_A,IDICHR_A,I_SO_A,ISFLIP_A,ICHKDIR_A,IPHI_A
     &,ITHETA_A,IE_A
        WRITE(IUO2,21) NPHI0,NTHETA0,NE,NPLAN,ISOM
      ELSEIF(JEL.EQ.2) THEN
        WRITE(IUO2,20) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        WRITE(IUO2,21) NPHI0,NTHETA0,NE,NPLAN,ISOM
      ENDIF
C
      DO JE=1,NE
        DO JANGLE=1,NANGLE
          IF(SPECTRO.EQ.'APC') THEN
            IF(JEL.EQ.1) THEN
              THETA=DTHETAA(JANGLE)
              PHI=DPHIA(JANGLE)
            ELSEIF(JEL.EQ.2) THEN
              THETA=DTHETA(JANGLE)
              PHI=DPHI(JANGLE)
            ENDIF
          ENDIF
C
          DO JPLAN=1,NPLAN
            IF(IDICHR.EQ.0) THEN
              IF((SPECTRO.NE.'APC').OR.(JEL.EQ.0)) THEN
                WRITE(IUO2,33) JPLAN,ECIN(JE),SUMR_1(JPLAN,JE,JANGLE),SU
     &MF_1(JPLAN,JE,JANGLE)
              ELSE
                WRITE(IUO2,34) JPLAN,THETA,PHI,ECIN(JE),SUMR_1(JPLAN,JE,
     &JANGLE),SUMF_1(JPLAN,JE,JANGLE)
              ENDIF
            ELSE
              IF((SPECTRO.NE.'APC').OR.(JEL.EQ.0)) THEN
                WRITE(IUO2,43) JPLAN,ECIN(JE),SUMR_1(JPLAN,JE,JANGLE),SU
     &MF_1(JPLAN,JE,JANGLE),SUMR_2(JPLAN,JE,JANGLE),SUMF_2(JPLAN,JE,JANG
     &LE)
              ELSE
                WRITE(IUO2,44) JPLAN,THETA,PHI,ECIN(JE),SUMR_1(JPLAN,JE,
     &JANGLE),SUMF_1(JPLAN,JE,JANGLE),SUMR_2(JPLAN,JE,JANGLE),SUMF_2(JPL
     &AN,JE,JANGLE)
              ENDIF
            ENDIF
          ENDDO
C
          IF(IDICHR.EQ.0) THEN
            IF((SPECTRO.NE.'APC').OR.(JEL.EQ.0)) THEN
              WRITE(IUO2,33) JVOL,ECIN(JE),SUMR_1(NPLAN+1,JE,JANGLE),SUM
     &F_1(NPLAN+1,JE,JANGLE)
              WRITE(IUO2,33) JTOT,ECIN(JE),SUMR_1(NPLAN+2,JE,JANGLE),SUM
     &F_1(NPLAN+2,JE,JANGLE)
            ELSE
              WRITE(IUO2,34) JVOL,THETA,PHI,ECIN(JE),SUMR_1(NPLAN+1,JE,J
     &ANGLE),SUMF_1(NPLAN+1,JE,JANGLE)
              WRITE(IUO2,34) JTOT,THETA,PHI,ECIN(JE),SUMR_1(NPLAN+2,JE,J
     &ANGLE),SUMF_1(NPLAN+2,JE,JANGLE)
            ENDIF
          ELSE
            IF((SPECTRO.NE.'APC').OR.(JEL.EQ.0)) THEN
              WRITE(IUO2,43) JVOL,ECIN(JE),SUMR_1(NPLAN+1,JE,JANGLE),SUM
     &F_1(NPLAN+1,JE,JANGLE),SUMR_2(NPLAN+1,JE,JANGLE),SUMF_2(NPLAN+1,JE
     &,JANGLE)
              WRITE(IUO2,43) JTOT,ECIN(JE),SUMR_1(NPLAN+2,JE,JANGLE),SUM
     &F_1(NPLAN+2,JE,JANGLE),SUMR_2(NPLAN+2,JE,JANGLE),SUMF_2(NPLAN+2,JE
     &,JANGLE)
            ELSE
              WRITE(IUO2,44) JVOL,THETA,PHI,ECIN(JE),SUMR_1(NPLAN+1,JE,J
     &ANGLE),SUMF_1(NPLAN+1,JE,JANGLE),SUMR_2(NPLAN+1,JE,JANGLE),SUMF_2(
     &NPLAN+1,JE,JANGLE)
              WRITE(IUO2,44) JTOT,THETA,PHI,ECIN(JE),SUMR_1(NPLAN+2,JE,J
     &ANGLE),SUMF_1(NPLAN+2,JE,JANGLE),SUMR_2(NPLAN+2,JE,JANGLE),SUMF_2(
     &NPLAN+2,JE,JANGLE)
            ENDIF
          ENDIF
C
        ENDDO
      ENDDO
C
   1  FORMAT(13X,I4)
   2  FORMAT(15X,F8.3,3X,F8.3,3X,E12.6)
   3  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
   4  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL ','IN 
     &THE WEIGHT_SUM SUBROUTINE - INCREASE NPM TO ',I3,'>>>>>>>>>>')
   5  FORMAT(6X,I1,1X,I3,3X,I3)
   8  FORMAT(I4,2X,I4,2X,I4,2X,I3,2X,I1)
   9  FORMAT(9(2X,I1),2X,I2)
  13  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,F6.2,2X,E12.6,2X,E
     &12.6)
  15  FORMAT(2X,A3,11X,A13)
  16  FORMAT(2X,A3,A5,1X,A3,2X,A13)
  18  FORMAT(I4,2X,I3,2X,I1)
  19  FORMAT(4(2X,I1))
  20   FORMAT(8(2X,I1))
  21   FORMAT(I4,2X,I4,2X,I4,2X,I3,2X,I1)
  23  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6,2X,E12.6,2X
     &,E12.6)
  24  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,F6.2,2X,E12.6,2X,E
     &12.6,2X,E12.6,2X,E12.6)
  33  FORMAT(2X,I3,2X,F8.2,2X,E12.6,2X,E12.6)
  34  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
  43  FORMAT(2X,I3,2X,F8.2,2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6)
  44  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6,2X,E12.6,2X
     &,E12.6)
C
      RETURN
C
      END
