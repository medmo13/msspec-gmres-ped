C
C=======================================================================
C
      SUBROUTINE STOP_TREAT(NFICHLEC,NPLAN,NEMET,NE,NTHETA,NTHETA_A,NPHI
     &,NPHI_A,ISOM,I_EXT,I_EXT_A,SPECTRO)
C
C   This subroutine stops the code before the long MS calculations
C     when the dimensioning NDIM_M of the treatment routines
C     (treat_aed,treat_apc,treat_phd,treat_xas) is insufficient.
C
C
C                                         Last modified : 06 Oct 2006
C
      USE DIM_MOD
      USE OUTUNITS_MOD
C
      CHARACTER*3 SPECTRO
C
C
      IF(ISOM.EQ.0) THEN
C
C   Photoelectron diffraction case
C
        IF(SPECTRO.EQ.'PHD') THEN
          IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
            NDP=NEMET*NTHETA*NPHI*NE
          ELSEIF(I_EXT.EQ.-1) THEN
            NDP=NEMET*NTHETA*NPHI*NE*2
          ELSEIF(I_EXT.EQ.2) THEN
            NDP=NEMET*NTHETA*NE
          ENDIF
          NTT=NPLAN*NDP
          IF(NTT.GT.NDIM_M) GOTO 10
          IF((NTHETA.GT.NTH_M).OR.(NPHI.GT.NPH_M)) GOTO 50
C
C   Auger electron diffraction case
C
        ELSEIF(SPECTRO.EQ.'AED') THEN
          IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
            NDP=NEMET*NTHETA_A*NPHI_A*NE
          ELSEIF(I_EXT_A.EQ.-1) THEN
            NDP=NEMET*NTHETA_A*NPHI_A*NE*2
          ELSEIF(I_EXT_A.EQ.2) THEN
            NDP=NEMET*NTHETA_A*NE
          ENDIF
          NTT=NPLAN*NDP
          IF(NTT.GT.NDIM_M) GOTO 20
          IF((NTHETA_A.GT.NTH_M).OR.(NPHI_A.GT.NPH_M)) GOTO 50
C
C   X-ray absorption case
C
        ELSEIF(SPECTRO.EQ.'XAS') THEN
          NDP=NEMET*NE
          NTT=NPLAN*NDP
          IF(NTT.GT.NDIM_M) GOTO 30
C
C   Auger Photoelectron coincidence spectroscopy case
C
        ELSEIF(SPECTRO.EQ.'APC') THEN
          IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
            IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
              NDP=NEMET*NTHETA*NPHI*NE*NTHETA_A*NPHI_A
            ELSEIF(I_EXT_A.EQ.-1) THEN
              NDP=NEMET*NTHETA*NPHI*NE*NTHETA_A*NPHI_A*2
            ELSEIF(I_EXT_A.EQ.2) THEN
              NDP=NEMET*NTHETA*NPHI*NE*NTHETA_A
            ENDIF
          ELSEIF(I_EXT.EQ.-1) THEN
            IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
              NDP=NEMET*NTHETA*NPHI*NE*2*NTHETA_A*NPHI_A
            ELSEIF(I_EXT_A.EQ.-1) THEN
              NDP=NEMET*NTHETA*NPHI*NE*2*NTHETA_A*NPHI_A*2
            ELSEIF(I_EXT_A.EQ.2) THEN
              NDP=NEMET*NTHETA*NPHI*NE*2*NTHETA_A
            ENDIF
          ELSEIF(I_EXT.EQ.2) THEN
            IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
              NDP=NEMET*NTHETA*NE*NTHETA_A*NPHI_A
            ELSEIF(I_EXT_A.EQ.-1) THEN
              NDP=NEMET*NTHETA*NE*NTHETA_A*NPHI_A*2
            ELSEIF(I_EXT_A.EQ.2) THEN
              NDP=NEMET*NTHETA*NE*NTHETA_A
            ENDIF
          ENDIF
          NTT=NPLAN*NDP
          IF(NTT.GT.NDIM_M) GOTO 40
          IF((NTHETA.GT.NTH_M).OR.(NPHI.GT.NPH_M)) GOTO 50
          IF((NTHETA_A.GT.NTH_M).OR.(NPHI_A.GT.NPH_M)) GOTO 50
          ENDIF
C
      ELSE
C
C   Photoelectron diffraction case
C
        IF(SPECTRO.EQ.'PHD') THEN
          IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
            NDP=NTHETA*NPHI*NE
          ELSEIF(I_EXT.EQ.-1) THEN
            NDP=NTHETA*NPHI*NE*2
          ELSEIF(I_EXT.EQ.2) THEN
            NDP=NTHETA*NE
          ENDIF
          NTT=NFICHLEC*NDP
          IF(NTT.GT.NDIM_M) GOTO 10
          IF((NTHETA.GT.NTH_M).OR.(NPHI.GT.NPH_M)) GOTO 50
C
C   Auger electron diffraction case
C
        ELSEIF(SPECTRO.EQ.'AED') THEN
          IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
            NDP=NTHETA_A*NPHI_A*NE
          ELSEIF(I_EXT_A.EQ.-1) THEN
            NDP=NTHETA_A*NPHI_A*NE*2
          ELSEIF(I_EXT_A.EQ.2) THEN
            NDP=NTHETA_A*NE
          ENDIF
          NTT=NFICHLEC*NDP
          IF(NTT.GT.NDIM_M) GOTO 20
          IF((NTHETA_A.GT.NTH_M).OR.(NPHI_A.GT.NPH_M)) GOTO 50
C
C   X-ray absorption case
C
        ELSEIF(SPECTRO.EQ.'XAS') THEN
          NDP=NE
          NTT=NFICHLEC*NDP
          IF(NTT.GT.NDIM_M) GOTO 30
C
C   Auger Photoelectron coincidence spectroscopy case
C
        ELSEIF(SPECTRO.EQ.'APC') THEN
          IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
            IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
              NDP=NTHETA*NPHI*NE*NTHETA_A*NPHI_A
            ELSEIF(I_EXT_A.EQ.-1) THEN
              NDP=NTHETA*NPHI*NE*NTHETA_A*NPHI_A*2
            ELSEIF(I_EXT_A.EQ.2) THEN
              NDP=NTHETA*NPHI*NE*NTHETA_A
            ENDIF
          ELSEIF(I_EXT.EQ.-1) THEN
            IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
              NDP=NTHETA*NPHI*NE*2*NTHETA_A*NPHI_A
            ELSEIF(I_EXT_A.EQ.-1) THEN
              NDP=NTHETA*NPHI*NE*2*NTHETA_A*NPHI_A*2
            ELSEIF(I_EXT_A.EQ.2) THEN
              NDP=NTHETA*NPHI*NE*2*NTHETA_A
            ENDIF
          ELSEIF(I_EXT.EQ.2) THEN
            IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
              NDP=NTHETA*NE*NTHETA_A*NPHI_A
            ELSEIF(I_EXT_A.EQ.-1) THEN
              NDP=NTHETA*NE*NTHETA_A*NPHI_A*2
            ELSEIF(I_EXT_A.EQ.2) THEN
              NDP=NTHETA*NE*NTHETA_A
            ENDIF
          ENDIF
          NTT=NFICHLEC*NDP
          IF(NTT.GT.NDIM_M) GOTO 40
          IF((NTHETA.GT.NTH_M).OR.(NPHI.GT.NPH_M)) GOTO 50
          IF((NTHETA_A.GT.NTH_M).OR.(NPHI_A.GT.NPH_M)) GOTO 50
        ENDIF
      ENDIF
C
      GOTO 5
C
  10  WRITE(IUO1,11) NTT
      STOP
  20  WRITE(IUO1,21) NTT
      STOP
  30  WRITE(IUO1,31) NTT
      STOP
  40  WRITE(IUO1,41) NTT
      STOP
  50  WRITE(IUO1,51) MAX(NTHETA,NPHI,NTHETA_A,NPHI_A)
      STOP
C
  11  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL',' IN 
     &THE    >>>>>>>>>>',/,8X,'<<<<<<<<<<  INCLUDE FILE ','FOR THE TREAT
     &_PHD SUBROUTINE   >>>>>>>>>>',/,8X,'<<<<<<<<<<          NDIM_M BE 
     &AT LEAST ',I8,'         >>>>>>>>>>')
  21  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL',' IN 
     &THE    >>>>>>>>>>',/,8X,'<<<<<<<<<<  INCLUDE FILE ','FOR THE TREAT
     &_AED SUBROUTINE   >>>>>>>>>>',/,8X,'<<<<<<<<<<          NDIM_M BE 
     &AT LEAST ',I8,'         >>>>>>>>>>')
  31  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL',' IN 
     &THE    >>>>>>>>>>',/,8X,'<<<<<<<<<<  INCLUDE FILE ','FOR THE TREAT
     &_XAS SUBROUTINE   >>>>>>>>>>',/,8X,'<<<<<<<<<<          NDIM_M BE 
     &AT LEAST ',I8,'         >>>>>>>>>>')
  41  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL',' IN 
     &THE    >>>>>>>>>>',/,8X,'<<<<<<<<<<  INCLUDE FILE ','FOR THE TREAT
     &_APC SUBROUTINE   >>>>>>>>>>',/,8X,'<<<<<<<<<<          NDIM_M BE 
     &AT LEAST ',I8,'         >>>>>>>>>>')
  51  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF NTH_M OR NPH_M TOO SMALL','
     &IN THE INCLUDE FILE - SHOULD BE AT LEAST ',I6,'  >>>>>>>>>>')
C
   5  RETURN
C
      END
