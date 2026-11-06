C
C=======================================================================
C
      SUBROUTINE TREAT_PHD(ISOM,NFICHLEC,JFICH,NP)
C
C  This routine sums up the calculations corresponding to different
C        absorbers or different planes when this has to be done
C             (parameter ISOM in the input data file).
C
C                                       Last modified : 24 Jan 2013
C
      USE DIM_MOD
      USE OUTUNITS_MOD
      USE TYPEXP_MOD , DUMMY => SPECTRO
      USE VALIN_MOD
      USE VALFIN_MOD
C
      PARAMETER(N_HEAD=5000,N_FILES=1000)
C
      CHARACTER*3 SPECTRO
C
      CHARACTER*13 OUTDATA
      CHARACTER*72 HEAD(N_HEAD,N_FILES)
C
      REAL TAB(NDIM_M,4)
      REAL ECIN(NE_M),DTHETA(NTH_M),DPHI(NPH_M)
C
C
      DATA JVOL,JTOT/0,-1/
C
      REWIND IUO2
C
C  Reading and storing the headers:
C
      NHEAD=0
      DO JLINE=1,N_HEAD
        READ(IUO2,888) HEAD(JLINE,JFICH)
        NHEAD=NHEAD+1
        IF(HEAD(JLINE,JFICH)(1:6).EQ.'      ') GOTO 333
      ENDDO
C
 333  CONTINUE
C
      READ(IUO2,15) SPECTRO,OUTDATA
      READ(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE,IPH_1
     &,I_EXT
C
      IF(I_EXT.EQ.2) THEN
         IPH_1=0
      ENDIF
C
      IF(ISOM.EQ.0) THEN
C
C........  ISOM = 0 : case of independent input files  .................
C
        READ(IUO2,1) NPLAN,NEMET,NTHETA,NPHI,NE
C
        IF(IPH_1.EQ.1) THEN
          N_FIXED=NPHI
          FIX0=PHI0
          FIX1=PHI1
          N_SCAN=NTHETA
        ELSE
          N_FIXED=NTHETA
          FIX0=THETA0
          FIX1=THETA1
          IF(STEREO.EQ.'YES') THEN
            NPHI=INT((PHI1-PHI0)*FLOAT(NTHETA-1)/(THETA1-THETA0)+0.0001)
     &+1
            IF(NTHETA*NPHI.GT.NPH_M) GOTO 37
          ENDIF
          N_SCAN=NPHI
        ENDIF
C
        IF(I_EXT.EQ.-1) THEN
          N_SCAN=2*N_SCAN
        ENDIF
C
        IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
          NDP=NEMET*NTHETA*NPHI*NE
        ELSEIF(I_EXT.EQ.-1) THEN
          NDP=NEMET*NTHETA*NPHI*NE*2
        ELSEIF(I_EXT.EQ.2) THEN
          NDP=NEMET*NTHETA*NE
          N_FIXED=NTHETA
          N_SCAN=NPHI
          IF((N_FIXED.GT.NTH_M).OR.(N_FIXED.GT.NPH_M)) GOTO 35
        ENDIF
C
        NTT=NPLAN*NDP
        IF(NTT.GT.NDIM_M) GOTO 5
C
        DO JPLAN=1,NPLAN
         DO JEMET=1,NEMET
          DO JE=1,NE
C
           DO J_FIXED=1,N_FIXED
            IF(N_FIXED.GT.1) THEN
              XINCRF=FLOAT(J_FIXED-1)*(FIX1-FIX0)/FLOAT(N_FIXED-1)
            ELSEIF(N_FIXED.EQ.1) THEN
              XINCRF=0.
            ENDIF
            IF(IPH_1.EQ.1) THEN
              JPHI=J_FIXED
            ELSE
              THETA=THETA0+XINCRF
              JTHETA=J_FIXED
              IF((ABS(THETA).GT.90.).AND.(I_EXT.NE.2)) GOTO 11
            ENDIF
            IF(STEREO.EQ.' NO') THEN
              N_SCAN_R=N_SCAN
            ELSE
              RTHETA=THETA*0.017453
              FIX_STEP=(FIX1-FIX0)/FLOAT(N_FIXED-1)
              N_SCAN_R=INT((PHI1-PHI0)*SIN(RTHETA)/FIX_STEP+0.0001)+1
            ENDIF
C
            DO J_SCAN=1,N_SCAN_R
             IF(IPH_1.EQ.1) THEN
              JTHETA=J_SCAN
             ELSE
              JPHI=J_SCAN
             ENDIF
C
             JLIN=(JPLAN-1)*NDP + (JEMET-1)*NE*N_FIXED*N_SCAN + (JE-1)*N
     &_FIXED*N_SCAN +(JTHETA-1)*NPHI + JPHI
C
             IF(I_EXT.LE.0) THEN
               IF(STEREO.EQ.' NO') THEN
                 JPHI2=JPHI
               ELSE
                 JPHI2=(JTHETA-1)*NPHI+JPHI
               ENDIF
             ELSE
               JPHI2=JTHETA
             ENDIF
C
             READ(IUO2,2) JPL
             IF(JPLAN.EQ.JPL) THEN
               BACKSPACE IUO2
               IF(IDICHR.EQ.0) THEN
                 READ(IUO2,2) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE
     &),TAB(JLIN,1),TAB(JLIN,2)
                 IF(I_EXT.EQ.-1) THEN
                   JLIN2=NTT+JLIN
                   READ(IUO2,25) TAB(JLIN2,1),TAB(JLIN2,2)
                 ENDIF
               ELSE
                 READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),ECIN(J
     &E),TAB(JLIN,1),TAB(JLIN,2),TAB(JLIN,3),TAB(JLIN,4)
                 IF(I_EXT.EQ.-1) THEN
                   JLIN2=NTT+JLIN
                   READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),ECIN
     &(JE),TAB(JLIN2,1),TAB(JLIN2,2),TAB(JLIN2,3),TAB(JLIN2,4)
                 ENDIF
               ENDIF
             ELSE
               BACKSPACE IUO2
               DO JL=JLIN,JPLAN*NDP
                 TAB(JL,1)=0.0
                 TAB(JL,2)=0.0
                 TAB(JL,3)=0.0
                 TAB(JL,4)=0.0
               ENDDO
               GOTO 10
             ENDIF
            ENDDO
           ENDDO
  11       CONTINUE
          ENDDO
         ENDDO
  10     CONTINUE
        ENDDO
C
        REWIND IUO2
C
C  Skipping the NHEAD lines of headers before rewriting:
C
        DO JLINE=1,NHEAD
          READ(IUO2,888) HEAD(JLINE,JFICH)
        ENDDO
C
        WRITE(IUO2,15) SPECTRO,OUTDATA
        WRITE(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        WRITE(IUO2,8) NPHI,NTHETA,NE,NPLAN,ISOM
C
        DO JE=1,NE
         DO JTHETA=1,NTHETA
          IF(STEREO.EQ.' NO') THEN
            NPHI_R=NPHI
          ELSE
            RTHETA=DTHETA(JTHETA)*0.017453
            FIX_STEP=(THETA1-THETA0)/FLOAT(NTHETA-1)
            NPHI_R=INT((PHI1-PHI0)*SIN(RTHETA)/FIX_STEP+0.0001)+1
            NPHI=INT((PHI1-PHI0)/FIX_STEP+0.0001)+1
          ENDIF
          DO JPHI=1,NPHI_R
           TOTDIF_1=0.
           TOTDIR_1=0.
           VOLDIF_1=0.
           VOLDIR_1=0.
           TOTDIF_2=0.
           TOTDIR_2=0.
           VOLDIF_2=0.
           VOLDIR_2=0.
           IF(I_EXT.EQ.-1) THEN
             TOTDIF2_1=0.
             TOTDIR2_1=0.
             VOLDIF2_1=0.
             VOLDIR2_1=0.
             TOTDIF2_2=0.
             TOTDIR2_2=0.
             VOLDIF2_2=0.
             VOLDIR2_2=0.
           ENDIF
C
           DO JPLAN=1,NPLAN
C
            SF_1=0.
            SR_1=0.
            SF_2=0.
            SR_2=0.
            IF(I_EXT.EQ.-1) THEN
              SF2_1=0.
              SR2_1=0.
              SF2_2=0.
              SR2_2=0.
            ENDIF
C
            DO JEMET=1,NEMET
             JLIN=(JPLAN-1)*NDP + (JEMET-1)*NE*NTHETA*NPHI + (JE-1)*NTHE
     &TA*NPHI +(JTHETA-1)*NPHI + JPHI
             SF_1=SF_1+TAB(JLIN,2)
             SR_1=SR_1+TAB(JLIN,1)
             IF(I_EXT.EQ.-1) THEN
               JLIN2=NTT+JLIN
               SF2_1=SF2_1+TAB(JLIN2,2)
               SR2_1=SR2_1+TAB(JLIN2,1)
             ENDIF
             IF(IDICHR.GE.1) THEN
               SF_2=SF_2+TAB(JLIN,4)
               SR_2=SR_2+TAB(JLIN,3)
               IF(I_EXT.EQ.-1) THEN
                 JLIN2=NTT+JLIN
                 SF2_2=SF2_2+TAB(JLIN2,4)
                 SR2_2=SR2_2+TAB(JLIN2,3)
               ENDIF
             ENDIF
            ENDDO
            IF(I_EXT.LE.0) THEN
               IF(STEREO.EQ.' NO') THEN
                 JPHI2=JPHI
               ELSE
                 JPHI2=(JTHETA-1)*NPHI+JPHI
               ENDIF
             ELSE
              JPHI2=JTHETA
            ENDIF
            IF(IDICHR.EQ.0) THEN
              WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),SR
     &_1,SF_1
              IF(I_EXT.EQ.-1) THEN
                WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),
     &SR2_1,SF2_1
              ENDIF
            ELSE
              WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),S
     &R_1,SF_1,SR_2,SF_2
              IF(I_EXT.EQ.-1) THEN
                WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE)
     &,SR2_1,SF2_1,SR2_2,SF2_2
              ENDIF
            ENDIF
            IF(JPLAN.GT.NONVOL(JFICH)) THEN
              VOLDIF_1=VOLDIF_1+SF_1
              VOLDIR_1=VOLDIR_1+SR_1
              IF(I_EXT.EQ.-1) THEN
                VOLDIF2_1=VOLDIF2_1+SF2_1
                VOLDIR2_1=VOLDIR2_1+SR2_1
              ENDIF
              IF(IDICHR.GE.1) THEN
                VOLDIF_2=VOLDIF_2+SF_2
                VOLDIR_2=VOLDIR_1+SR_2
                IF(I_EXT.EQ.-1) THEN
                  VOLDIF2_2=VOLDIF2_2+SF2_2
                  VOLDIR2_2=VOLDIR2_1+SR2_2
                ENDIF
              ENDIF
            ENDIF
            TOTDIF_1=TOTDIF_1+SF_1
            TOTDIR_1=TOTDIR_1+SR_1
            IF(I_EXT.EQ.-1) THEN
              TOTDIF2_1=TOTDIF2_1+SF2_1
              TOTDIR2_1=TOTDIR2_1+SR2_1
            ENDIF
            IF(IDICHR.GE.1) THEN
              TOTDIF_2=TOTDIF_2+SF_2
              TOTDIR_2=TOTDIR_2+SR_2
              IF(I_EXT.EQ.-1) THEN
                TOTDIF2_2=TOTDIF2_2+SF2_2
                TOTDIR2_2=TOTDIR2_2+SR2_2
              ENDIF
            ENDIF
           ENDDO
           IF(IDICHR.EQ.0) THEN
             WRITE(IUO2,3) JVOL,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),VOLD
     &IR_1,VOLDIF_1
             IF(I_EXT.EQ.-1) THEN
               WRITE(IUO2,3) JVOL,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),VO
     &LDIR2_1,VOLDIF2_1
             ENDIF
             WRITE(IUO2,3) JTOT,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),TOTD
     &IR_1,TOTDIF_1
             IF(I_EXT.EQ.-1) THEN
               WRITE(IUO2,3) JTOT,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),TO
     &TDIR2_1,TOTDIF2_1
             ENDIF
           ELSE
             WRITE(IUO2,23) JVOL,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),VOL
     &DIR_1,VOLDIF_1,VOLDIR_2,VOLDIF_2
            IF(I_EXT.EQ.-1) THEN
               WRITE(IUO2,23) JVOL,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),V
     &OLDIR2_1,VOLDIF2_1,VOLDIR2_2,VOLDIF2_2
             ENDIF
             WRITE(IUO2,23) JTOT,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),TOT
     &DIR_1,TOTDIF_1,TOTDIR_2,TOTDIF_2
             IF(I_EXT.EQ.-1) THEN
               WRITE(IUO2,23) JTOT,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),T
     &OTDIR2_1,TOTDIF2_1,TOTDIR2_2,TOTDIF2_2
             ENDIF
           ENDIF
          ENDDO
         ENDDO
        ENDDO
C
      ELSE
C
C........ ISOM not= 0 : multiple input files to be summed up  ..........
C
        READ(IUO2,7) NTHETA,NPHI,NE
C
        IF(IPH_1.EQ.1) THEN
          N_FIXED=NPHI
          FIX0=PHI0
          FIX1=PHI1
          N_SCAN=NTHETA
        ELSE
          N_FIXED=NTHETA
          FIX0=THETA0
          FIX1=THETA1
          IF(STEREO.EQ.'YES') THEN
            NPHI=INT((PHI1-PHI0)*FLOAT(NTHETA-1)/(THETA1-THETA0)+0.0001)
     &+1
            IF(NTHETA*NPHI.GT.NPH_M) GOTO 37
          ENDIF
          N_SCAN=NPHI
        ENDIF
C
        IF(I_EXT.EQ.-1) THEN
          N_SCAN=2*N_SCAN
        ENDIF
C
        IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
          NDP=NTHETA*NPHI*NE
        ELSEIF(I_EXT.EQ.-1) THEN
          NDP=NTHETA*NPHI*NE*2
        ELSEIF(I_EXT.EQ.2) THEN
          NDP=NTHETA*NE
          N_FIXED=NTHETA
          N_SCAN=NPHI
          IF((N_FIXED.GT.NTH_M).OR.(N_FIXED.GT.NPH_M)) GOTO 35
        ENDIF
C
        NTT=NFICHLEC*NDP
        IF(NTT.GT.NDIM_M) GOTO 5
C
        IF(ISOM.EQ.1) THEN
          NPLAN=NP
          NF=NP
        ELSEIF(ISOM.EQ.2) THEN
          NEMET=NFICHLEC
          NF=NFICHLEC
          NPLAN=1
        ENDIF
C
        DO JF=1,NF
C
C  Reading the headers for each file:
C
          IF(JF.GT.1) THEN
            DO JLINE=1,NHEAD
              READ(IUO2,888) HEAD(JLINE,JF)
            ENDDO
          ENDIF
C
          DO JE=1,NE
C
           DO J_FIXED=1,N_FIXED
            IF(N_FIXED.GT.1) THEN
              XINCRF=FLOAT(J_FIXED-1)*(FIX1-FIX0)/FLOAT(N_FIXED-1)
            ELSEIF(N_FIXED.EQ.1) THEN
              XINCRF=0.
            ENDIF
            IF(IPH_1.EQ.1) THEN
              JPHI=J_FIXED
            ELSE
              THETA=THETA0+XINCRF
              JTHETA=J_FIXED
              IF((ABS(THETA).GT.90.).AND.(I_EXT.NE.2)) GOTO 12
            ENDIF
            IF(STEREO.EQ.' NO') THEN
              N_SCAN_R=N_SCAN
            ELSE
              RTHETA=THETA*0.017453
              FIX_STEP=(FIX1-FIX0)/FLOAT(N_FIXED-1)
              N_SCAN_R=INT((PHI1-PHI0)*SIN(RTHETA)/FIX_STEP+0.0001)+1
            ENDIF
C
            DO J_SCAN=1,N_SCAN_R
             IF(IPH_1.EQ.1) THEN
              JTHETA=J_SCAN
             ELSE
              JPHI=J_SCAN
             ENDIF
C
             JLIN=(JF-1)*NDP + (JE-1)*N_FIXED*N_SCAN +(JTHETA-1)*NPHI + 
     &JPHI
             IF(I_EXT.LE.0) THEN
               IF(STEREO.EQ.' NO') THEN
                 JPHI2=JPHI
               ELSE
                 JPHI2=(JTHETA-1)*NPHI+JPHI
               ENDIF
             ELSE
               JPHI2=JTHETA
             ENDIF
C
             IF(ISOM.EQ.1) THEN
               READ(IUO2,2) JPL
               IF(JF.EQ.JPL) THEN
                 BACKSPACE IUO2
                 IF(IDICHR.EQ.0) THEN
                   READ(IUO2,2) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),ECIN(
     &JE),TAB(JLIN,1),TAB(JLIN,2)
                   IF(I_EXT.EQ.-1) THEN
                     JLIN2=NTT+JLIN
                     READ(IUO2,25) TAB(JLIN2,1),TAB(JLIN2,2)
                   ENDIF
                 ELSE
                   READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),ECIN
     &(JE),TAB(JLIN,1),TAB(JLIN,2),TAB(JLIN,3),TAB(JLIN,4)
                   IF(I_EXT.EQ.-1) THEN
                     JLIN2=NTT+JLIN
                     READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),EC
     &IN(JE),TAB(JLIN2,1),TAB(JLIN2,2),TAB(JLIN2,3),TAB(JLIN2,4)
                   ENDIF
                 ENDIF
               ELSE
                 BACKSPACE IUO2
                 DO JLINE=1,NHEAD
                   BACKSPACE IUO2
                 ENDDO
                 DO JL=JLIN,JF*NDP
                   TAB(JL,1)=0.0
                   TAB(JL,2)=0.0
                   TAB(JL,3)=0.0
                   TAB(JL,4)=0.0
                 ENDDO
                 GOTO 13
               ENDIF
             ELSEIF(ISOM.EQ.2) THEN
               IF(IDICHR.EQ.0) THEN
                 READ(IUO2,2) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE
     &),TAB(JLIN,1),TAB(JLIN,2)
                 IF(I_EXT.EQ.-1) THEN
                   JLIN2=NTT+JLIN
                   READ(IUO2,25) TAB(JLIN2,1),TAB(JLIN2,2)
                 ENDIF
               ELSE
                 READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),ECIN(J
     &E),TAB(JLIN,1),TAB(JLIN,2),TAB(JLIN,3),TAB(JLIN,4)
                 IF(I_EXT.EQ.-1) THEN
                   JLIN2=NTT+JLIN
                   READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),ECIN
     &(JE),TAB(JLIN2,1),TAB(JLIN2,2),TAB(JLIN2,3),TAB(JLIN2,4)
                   ENDIF
               ENDIF
             ENDIF
            ENDDO
 12         CONTINUE
           ENDDO
          ENDDO
 13       CONTINUE
        ENDDO
C
        REWIND IUO2
C
C  Writing the headers:
C
        DO JLINE=1,2
          WRITE(IUO2,888) HEAD(JLINE,1)
        ENDDO
        DO JF=1,NFICHLEC
          DO JLINE=3,6
            WRITE(IUO2,888) HEAD(JLINE,JF)
          ENDDO
          WRITE(IUO2,888) HEAD(2,JF)
        ENDDO
        DO JLINE=7,NHEAD
          WRITE(IUO2,888) HEAD(JLINE,1)
        ENDDO
C
        WRITE(IUO2,15) SPECTRO,OUTDATA
        WRITE(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        WRITE(IUO2,8) NPHI,NTHETA,NE,NPLAN,ISOM
C
        IF(ISOM.EQ.1) THEN
C
          DO JE=1,NE
C
           DO JTHETA=1,NTHETA
            IF(STEREO.EQ.' NO') THEN
              NPHI_R=NPHI
            ELSE
              RTHETA=DTHETA(JTHETA)*0.017453
              FIX_STEP=(THETA1-THETA0)/FLOAT(NTHETA-1)
              NPHI_R=INT((PHI1-PHI0)*SIN(RTHETA)/FIX_STEP+0.0001)+1
              NPHI=INT((PHI1-PHI0)/FIX_STEP+0.0001)+1
            ENDIF
            DO JPHI=1,NPHI_R
C
             TOTDIF_1=0.
             TOTDIR_1=0.
             VOLDIF_1=0.
             VOLDIR_1=0.
             TOTDIF_2=0.
             TOTDIR_2=0.
             VOLDIF_2=0.
             VOLDIR_2=0.
             IF(I_EXT.EQ.-1) THEN
               TOTDIF2_1=0.
               TOTDIR2_1=0.
               VOLDIF2_1=0.
               VOLDIR2_1=0.
               TOTDIF2_2=0.
               TOTDIR2_2=0.
               VOLDIF2_2=0.
               VOLDIR2_2=0.
             ENDIF
C
             DO JPLAN=1,NPLAN
              JF=JPLAN
C
              JLIN=(JF-1)*NDP + (JE-1)*NTHETA*NPHI +(JTHETA-1)*NPHI + JP
     &HI
C
              SR_1=TAB(JLIN,1)
              SF_1=TAB(JLIN,2)
              IF(I_EXT.EQ.-1) THEN
                 JLIN2=NTT+JLIN
                 SF2_1=TAB(JLIN2,2)
                 SR2_1=TAB(JLIN2,1)
              ENDIF
              IF(I_EXT.LE.0) THEN
                IF(STEREO.EQ.' NO') THEN
                  JPHI2=JPHI
                ELSE
                  JPHI2=(JTHETA-1)*NPHI+JPHI
                ENDIF
              ELSE
                JPHI2=JTHETA
              ENDIF
              IF(IDICHR.EQ.0) THEN
                WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),
     &SR_1,SF_1
                IF(I_EXT.EQ.-1) THEN
                  WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE
     &),SR2_1,SF2_1
                ENDIF
              ELSE
                SR_2=TAB(JLIN,3)
                SF_2=TAB(JLIN,4)
                IF(I_EXT.EQ.-1) THEN
                  JLIN2=NTT+JLIN
                  SF2_2=TAB(JLIN2,4)
                  SR2_2=TAB(JLIN2,3)
                ENDIF
                WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE)
     &,SR_1,SF_1,SR_2,SF_2
                IF(I_EXT.EQ.-1) THEN
                  WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),ECIN(J
     &E),SR2_1,SF2_1,SR2_2,SF2_2
                ENDIF
              ENDIF
              IF(NONVOL(JPLAN).EQ.0) THEN
                VOLDIF_1=VOLDIF_1+SF_1
                VOLDIR_1=VOLDIR_1+SR_1
                IF(I_EXT.EQ.-1) THEN
                  VOLDIF2_1=VOLDIF2_1+SF2_1
                  VOLDIR2_1=VOLDIR2_1+SR2_1
                ENDIF
                IF(IDICHR.GE.1) THEN
                  VOLDIF_2=VOLDIF_2+SF_2
                  VOLDIR_2=VOLDIR_2+SR_2
                  IF(I_EXT.EQ.-1) THEN
                    VOLDIF2_2=VOLDIF2_2+SF2_2
                    VOLDIR2_2=VOLDIR2_1+SR2_2
                  ENDIF
                ENDIF
              ENDIF
              TOTDIF_1=TOTDIF_1+SF_1
              TOTDIR_1=TOTDIR_1+SR_1
              IF(I_EXT.EQ.-1) THEN
                TOTDIF2_1=TOTDIF2_1+SF2_1
                TOTDIR2_1=TOTDIR2_1+SR2_1
              ENDIF
              IF(IDICHR.GE.1) THEN
                TOTDIF_2=TOTDIF_2+SF_2
                TOTDIR_2=TOTDIR_2+SR_2
                IF(I_EXT.EQ.-1) THEN
                  TOTDIF2_2=TOTDIF2_2+SF2_2
                  TOTDIR2_2=TOTDIR2_2+SR2_2
                ENDIF
              ENDIF
             ENDDO
C
             IF(IDICHR.EQ.0) THEN
               WRITE(IUO2,3) JVOL,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),VO
     &LDIR_1,VOLDIF_1
               IF(I_EXT.EQ.-1) THEN
                 WRITE(IUO2,3) JVOL,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),
     &VOLDIR2_1,VOLDIF2_1
               ENDIF
               WRITE(IUO2,3) JTOT,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),TO
     &TDIR_1,TOTDIF_1
               IF(I_EXT.EQ.-1) THEN
                 WRITE(IUO2,3) JTOT,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),
     &TOTDIR2_1,TOTDIF2_1
               ENDIF
             ELSE
               WRITE(IUO2,23) JVOL,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),V
     &OLDIR_1,VOLDIF_1,VOLDIR_2,VOLDIF_2
               IF(I_EXT.EQ.-1) THEN
                 WRITE(IUO2,23) JVOL,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE)
     &,VOLDIR2_1,VOLDIF2_1,VOLDIR2_2,VOLDIF2_2
               ENDIF
               WRITE(IUO2,23) JTOT,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),T
     &OTDIR_1,TOTDIF_1,TOTDIR_2,TOTDIF_2
               IF(I_EXT.EQ.-1) THEN
                 WRITE(IUO2,23) JTOT,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE)
     &,TOTDIR2_1,TOTDIF2_1,TOTDIR2_2,TOTDIF2_2
               ENDIF
             ENDIF
C
            ENDDO
           ENDDO
          ENDDO
        ELSEIF(ISOM.EQ.2) THEN
          DO JE=1,NE
C
           DO JTHETA=1,NTHETA
            IF(STEREO.EQ.' NO') THEN
              NPHI_R=NPHI
            ELSE
              RTHETA=DTHETA(JTHETA)*0.017453
              FIX_STEP=(THETA1-THETA0)/FLOAT(NTHETA-1)
              NPHI_R=INT((PHI1-PHI0)*SIN(RTHETA)/FIX_STEP+0.0001)+1
              NPHI=INT((PHI1-PHI0)/FIX_STEP+0.0001)+1
            ENDIF
            DO JPHI=1,NPHI_R
C
              SF_1=0.
              SR_1=0.
              SF_2=0.
              SR_2=0.
              IF(I_EXT.EQ.-1) THEN
                SF2_1=0.
                SR2_1=0.
                SF2_2=0.
                SR2_2=0.
              ENDIF
C
              DO JEMET=1,NEMET
               JF=JEMET
C
               JLIN=(JF-1)*NDP + (JE-1)*NTHETA*NPHI +(JTHETA-1)*NPHI + J
     &PHI
C
               SF_1=SF_1+TAB(JLIN,2)
               SR_1=SR_1+TAB(JLIN,1)
               IF(I_EXT.EQ.-1) THEN
                 JLIN2=NTT+JLIN
                 SF2_1=SF2_1+TAB(JLIN2,2)
                 SR2_1=SR2_1+TAB(JLIN2,1)
               ENDIF
               IF(IDICHR.GE.1) THEN
                 SF_2=SF_2+TAB(JLIN,4)
                 SR_2=SR_2+TAB(JLIN,3)
                 IF(I_EXT.EQ.-1) THEN
                   JLIN2=NTT+JLIN
                   SF2_2=SF2_2+TAB(JLIN2,4)
                   SR2_2=SR2_2+TAB(JLIN2,3)
                 ENDIF
               ENDIF
              ENDDO
              IF(I_EXT.LE.0) THEN
                IF(STEREO.EQ.' NO') THEN
                  JPHI2=JPHI
                ELSE
                  JPHI2=(JTHETA-1)*NPHI+JPHI
                ENDIF
              ELSE
                JPHI2=JTHETA
              ENDIF
              IF(IDICHR.EQ.0) THEN
                WRITE(IUO2,3) JPL,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),SR
     &_1,SF_1
                IF(I_EXT.EQ.-1) THEN
                  WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE
     &),SR2_1,SF2_1
                ENDIF
              ELSE
                WRITE(IUO2,23) JPL,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),S
     &R_1,SF_1,SR_2,SF_2
                IF(I_EXT.EQ.-1) THEN
                  WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),ECIN(J
     &E),SR2_1,SF2_1,SR2_2,SF2_2
                ENDIF
              ENDIF
            ENDDO
           ENDDO
          ENDDO
        ENDIF
      ENDIF
C
      GOTO 6
C
   5  WRITE(IUO1,4)
      STOP
  35  WRITE(IUO1,36) N_FIXED
      STOP
  37  WRITE(IUO1,38) NTHETA*NPHI
      STOP
C
   1  FORMAT(2X,I3,2X,I2,2X,I4,2X,I4,2X,I4)
   2  FORMAT(2X,I3,2X,I2,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
   3  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
   4  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL ','IN 
     &THE TREAT_PHD SUBROUTINE - INCREASE NDIM_M ','>>>>>>>>>>')
   7  FORMAT(I4,2X,I4,2X,I4)
   8  FORMAT(I4,2X,I4,2X,I4,2X,I3,2X,I1)
   9  FORMAT(9(2X,I1),2X,I2)
  15  FORMAT(2X,A3,11X,A13)
  22  FORMAT(2X,I3,2X,I2,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6,2X,E1
     &2.6,2X,E12.6)
  23  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6,2X,E12.6,2X
     &,E12.6)
  25  FORMAT(37X,E12.6,2X,E12.6)
  36  FORMAT(//,4X,'<<<<<<<<<<  DIMENSION OF NTH_M OR NPH_M TOO SMALL ',
     &'IN THE INCLUDE FILE >>>>>>>>>>',/,4X,'<<<<<<<<<<                 
     &SHOULD BE AT LEAST ',I6,'                 >>>>>>>>>>')
  38  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF NPH_M TOO SMALL ','IN THE I
     &NCLUDE FILE >>>>>>>>>>',/,8X,'<<<<<<<<<<             SHOULD BE AT 
     &LEAST ',I6,'             >>>>>>>>>>')
 888  FORMAT(A72)
C
   6  RETURN
C
      END
