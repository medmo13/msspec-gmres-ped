C
C=======================================================================
C
      SUBROUTINE READ_DATA(ICOM,NFICHLEC,JFICH,ITRTL,*,*,*,*,*,*,*,*,*,*
     &,*,*,*)
C
C  This subroutine reads the input data from unit ICOM and writes
C     them in the control file IUO1. Then, it stores the data in
C                     the various COMMON blocks
C
C                                        Last modified : 26 Apr 2013
C
      USE DIM_MOD
C
      USE ADSORB_MOD
      USE AMPLI_MOD
      USE APPROX_MOD
      USE ATOMS_MOD
      USE AUGER_MOD
      USE BASES_MOD
      USE COEFRLM_MOD
      USE CONVACC_MOD
      USE CONVTYP_MOD
      USE C_G_MOD
      USE C_G_A_MOD
      USE C_G_M_MOD
      USE CRANGL_MOD
      USE DEBWAL_MOD , T => TEMP
      USE DEXPFAC2_MOD
      USE DFACTSQ_MOD
      USE EIGEN_MOD
      USE EXAFS_MOD
      USE EXPFAC_MOD
      USE EXPFAC2_MOD
      USE EXPROT_MOD
      USE FACTSQ_MOD
      USE FDIF_MOD
      USE FIXSCAN_MOD
      USE FIXSCAN_A_MOD
      USE HEADER_MOD , AUGER1 => AUGER
      USE INDAT_MOD
      USE INFILES_MOD
      USE INUNITS_MOD
      USE INIT_A_MOD
      USE INIT_J_MOD
      USE INIT_L_MOD
      USE INIT_M_MOD
      USE LIMAMA_MOD
      USE LINLBD_MOD
      USE LOGAMAD_MOD
      USE LPMOY_MOD , XM => XMTA, RH => RHOTA
      USE MILLER_MOD
      USE MOYEN_MOD
      USE MOYEN_A_MOD
      USE OUTFILES_MOD
      USE OUTUNITS_MOD
      USE PARCAL_MOD
      USE PARCAL_A_MOD
      USE RA_MOD
      USE RELADS_MOD
      USE RELAX_MOD
      USE RENORM_MOD
      USE RESEAU_MOD
      USE SPECTRUM_MOD
      USE SPIN_MOD
      USE TESTS_MOD
      USE TYPCAL_MOD
      USE TYPCAL_A_MOD
      USE TYPEM_MOD
      USE TYPEXP_MOD
      USE VALIN_MOD
      USE VALIN_AV_MOD
      USE VALFIN_MOD
      USE VALEX_A_MOD
      USE XMRHO_MOD
C
C
C
C
      REAL*8 J1,J2,MJ1,MJ2,MJ3,JJ,DXDEN,DEXPF
      REAL*8 JJ_MIN,JJ_MAX,JJ12,JL12,SMALL,SQPI
C
      REAL TEXTE1(10),TEXTE2(10),TEXTE3(10)
      REAL TEXTE4(10),TEXTE5(10),TEXTE6(10)
      REAL TEXTE6B(10),TEXTE7(10)
      REAL THFWD(NATP_M),THBWD(NATP_M),GLG(0:N_GAUNT),NJ(0:N_GAUNT)
      REAL ALPHAR,BETAR,RACC
C
C
C
      DOUBLE PRECISION FACT1L,FACT2L
C
C
C
C
C
      CHARACTER*7 TESLEC,RIEN
C
C
      CHARACTER*3 CODRES(8),CODCTR(7),CRIST,CENTR,UNLENGTH
C
C
      CHARACTER*1 EDGE_C,EDGE_I,EDGE_A,MULT
      DATA CODRES/'CUB','TET','ORB','MNC','TCN','TRG','HEX','EXT'/
      DATA CODCTR/'P','I','F','R','A','B','C'/
      DATA PIS180,BOHR/0.017453,0.529177/
      DATA SQPI,SMALL /1.772453850906D0,1.D-6/
C
      I_EXT=0
      I_EXT_A=0
      IVG0=0
      IRET=0
      NCRIST=0
      NCENTR=0
      I_SO=0
      DO I=1,10
        PCREL(I)=0.
      ENDDO
      STEREO=' NO'
C
C
C..........   Reading of the input data in unit ICOM   ..........
C
C
      READ(ICOM,1) RIEN
      READ(ICOM,2) TEXTE1
      READ(ICOM,1) RIEN
      READ(ICOM,1) RIEN
      READ(ICOM,2) TEXTE2
      READ(ICOM,1) RIEN
C
      READ(ICOM,3) CRIST,CENTR,IBAS,NAT
      READ(ICOM,4) A,BSURA,CSURA,UNIT
C
      IF(IBAS.EQ.0) THEN
        DO JLINE=1,100
          READ(ICOM,5) TESLEC
          IF(TESLEC.EQ.'SPECTRO') THEN
            BACKSPACE ICOM
            BACKSPACE ICOM
            BACKSPACE ICOM
            GOTO 600
          ENDIF
        ENDDO
      ENDIF
C
      READ(ICOM,6) ALPHAD,BETAD,GAMMAD
      READ(ICOM,7) IH,IK,II,IL
      READ(ICOM,8) NIV,COUPUR,ITEST,IESURF
      IF(NAT.GT.1) THEN
        DO I=1,NAT
          J=3*(I-1)
          READ(ICOM,9) ATBAS(1+J),ATBAS(2+J),ATBAS(3+J),CHEM(I),NZAT(I)
        ENDDO
      ELSE
        READ(ICOM,9) X1,Y1,Z1,CHEM(1),NZA
      ENDIF
C
      READ(ICOM,5) TESLEC
      IF(TESLEC.EQ.'VECBAS ') THEN
        BACKSPACE ICOM
      ELSE
        IRET=10
        GOTO 605
      ENDIF
C
      DO I=1,8
        IF(CRIST.EQ.CODRES(I)) NCRIST=I
        IF(I.NE.8) THEN
          IF(CENTR.EQ.CODCTR(I)) NCENTR=I
        ENDIF
      ENDDO
      IF((NCRIST.EQ.0).OR.(NCENTR.EQ.0)) THEN
        IRET=1
        GOTO 605
      ENDIF
C
      IF(NCRIST.EQ.8) THEN
        DO I=1,3
          J=3*(I-1)
          IVN(I)=1
          READ(ICOM,9) VECBAS(1+J),VECBAS(2+J),VECBAS(3+J)
          IF(ABS(VECBAS(1+J)).LT.0.0001) THEN
            IF(ABS(VECBAS(2+J)).LT.0.0001) THEN
              IF(ABS(VECBAS(3+J)).LT.0.0001) THEN
                IVG0=IVG0+1
                IVN(I)=0
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ELSE
        READ(ICOM,9) X3,Y3,Z3
        READ(ICOM,9) X4,Y4,Z4
        READ(ICOM,9) X5,Y5,Z5
      ENDIF
      READ(ICOM,10) IREL,NREL,(PCREL(I),I=1,2)
      IF(IREL.EQ.1) THEN
        IF(NREL.GT.2) THEN
          NLIGNE=INT(FLOAT(NREL-2)/4.)+1
          DO J=1,NLIGNE
            READ(ICOM,11) (PCREL(I),I=1,4)
          ENDDO
        ENDIF
        IF(NREL.GT.10) THEN
          IRET=4
          GOTO 605
        ENDIF
      ELSEIF(IREL.EQ.0) THEN
        NREL=0
      ENDIF
      IF(NREL.EQ.0) THEN
        DO JREL=1,10
          PCREL(JREL)=0.
        ENDDO
      ENDIF
      READ(ICOM,12) OMEGAD1,OMEGAD2,IADS
C
      READ(ICOM,1) RIEN
 600  READ(ICOM,2) TEXTE3
      READ(ICOM,1) RIEN
C
      READ(ICOM,13) SPECTRO,ISPIN,IDICHR,IPOL
      READ(ICOM,44) I_AMP
C
      IF(SPECTRO.EQ.'PHD') THEN
        INTERACT='DIPOLAR'
      ELSEIF(SPECTRO.EQ.'LED') THEN
        INTERACT='NOINTER'
      ELSEIF(SPECTRO.EQ.'XAS') THEN
        INTERACT='DIPOLAR'
      ELSEIF(SPECTRO.EQ.'AED') THEN
        INTERACT='COULOMB'
      ELSEIF(SPECTRO.EQ.'APC') THEN
        INTERACT='DIPCOUL'
      ELSEIF(SPECTRO.EQ.'EIG') THEN
        INTERACT='DIPOLAR'
      ENDIF
C
      IF((IPOL.EQ.0).AND.(IDICHR.GT.0)) THEN
        PRINT 513
        STOP
      ENDIF
      IF((IDICHR.EQ.2).AND.(ISPIN.EQ.0)) THEN
        PRINT 514
        STOP
      ENDIF
C
      IF(ISPIN.EQ.0) THEN
        NSPIN2=1
        NSPIN=1
      ELSEIF(ISPIN.EQ.1) THEN
        NSPIN2=4
        NSPIN=2
      ENDIF
C
      IF(SPECTRO.EQ.'LED') THEN
        DO JLINE=1,10
          READ(ICOM,1) RIEN
        ENDDO
        GOTO 607
      ELSEIF(SPECTRO.EQ.'XAS') THEN
        IF(IDICHR.GT.1) THEN
          PRINT 512
          STOP
        ENDIF
        DO JLINE=1,19
          READ(ICOM,1) RIEN
        ENDDO
        GOTO 602
      ELSEIF(SPECTRO.EQ.'AED') THEN
        DO JLINE=1,24
          READ(ICOM,1) RIEN
        ENDDO
        GOTO 603
      ELSEIF(SPECTRO.EQ.'EIG') THEN
        DO JLINE=1,34
          READ(ICOM,1) RIEN
        ENDDO
        GOTO 608
      ENDIF
C
      IF((SPECTRO.EQ.'PHD').OR.(SPECTRO.EQ.'APC')) THEN
        READ(ICOM,1) RIEN
        READ(ICOM,2) TEXTE4
        READ(ICOM,1) RIEN
C
        READ(ICOM,20) NI,NLI,S_O,INITL,I_SO
C
        IF((NLI.EQ.'s').OR.(NLI.EQ.'S')) THEN
          LI=0
        ELSEIF((NLI.EQ.'p').OR.(NLI.EQ.'P')) THEN
          LI=1
        ELSEIF((NLI.EQ.'d').OR.(NLI.EQ.'D')) THEN
          LI=2
        ELSEIF((NLI.EQ.'f').OR.(NLI.EQ.'F')) THEN
          LI=3
        ELSEIF((NLI.EQ.'g').OR.(NLI.EQ.'G')) THEN
          LI=4
        ELSE
          IRET=5
          GOTO 605
        ENDIF
        IF(LI.GT.LI_M) THEN
          IRET=6
          GOTO 605
        ENDIF
        IF(I_SO.EQ.0) THEN
          S_O='   '
        ELSEIF(I_SO.EQ.1) THEN
          IF(S_O.EQ.'1/2') THEN
            IF(LI.GT.1) IRET=7
          ELSEIF(S_O.EQ.'3/2') THEN
            IF((LI.LT.1).OR.(LI.GT.2)) IRET=7
          ELSEIF(S_O.EQ.'5/2') THEN
            IF((LI.LT.2).OR.(LI.GT.3)) IRET=7
          ELSEIF(S_O.EQ.'7/2') THEN
            IF((LI.LT.3).OR.(LI.GT.4)) IRET=7
          ELSEIF(S_O.EQ.'9/2') THEN
            IF(LI.NE.4) IRET=7
          ENDIF
        ELSEIF(I_SO.EQ.2) THEN
          S_O='   '
        ENDIF
C
        READ(ICOM,14) IPHI,ITHETA,IE,IFTHET
        READ(ICOM,15) NPHI,NTHETA,NE,NFTHET
        READ(ICOM,16) PHI0,THETA0,E0,R1
        READ(ICOM,16) PHI1,THETA1,EFIN,R2
        READ(ICOM,17) THLUM,PHILUM,ELUM
        READ(ICOM,18) IMOD,IMOY,ACCEPT,ICHKDIR
C
        DO JLINE=1,9
          READ(ICOM,1) RIEN
        ENDDO
      ENDIF
C
 607  IF(SPECTRO.EQ.'LED') THEN
        READ(ICOM,1) RIEN
        READ(ICOM,2) TEXTE4
        READ(ICOM,1) RIEN
C
        READ(ICOM,14) IPHI,ITHETA,IE,IFTHET
        READ(ICOM,15) NPHI,NTHETA,NE,NFTHET
        READ(ICOM,16) PHI0,THETA0,E0,R1
        READ(ICOM,16) PHI1,THETA1,EFIN,R2
        READ(ICOM,17) TH_INI,PHI_INI
        READ(ICOM,18) IMOD,IMOY,ACCEPT,ICHKDIR
C
        THLUM=TH_INI
        PHILUM=PHI_INI
        ELUM=0.
        IDICHR=0
        INITL=0
      ENDIF
C
      IF(SPECTRO.NE.'XAS') THEN
        IF(IPHI.EQ.-1) THEN
          IPHI=1
          I_EXT=0
          ICHKDIR=0
          STEREO='YES'
          IF(ABS(PHI1-PHI0).LT.0.0001) THEN
            PHI0=0.
            PHI1=360.
            NPHI=361
          ENDIF
          IF(ABS(THETA1-THETA0).LT.0.0001) THEN
            THETA0=0.
            THETA1=88.
            NTHETA=89
          ENDIF
        ELSEIF(IPHI.EQ.2) THEN
          IPHI=1
          I_EXT=1
        ELSEIF(IPHI.EQ.3) THEN
          IPHI=1
          I_EXT=-1
        ELSEIF(ITHETA.EQ.2) THEN
          ITHETA=1
          I_EXT=1
        ELSEIF(ITHETA.EQ.3) THEN
          ITHETA=1
          I_EXT=-1
        ELSEIF(IE.EQ.2) THEN
          IE=1
          I_EXT=1
        ELSEIF(IE.EQ.3) THEN
          IE=1
          I_EXT=-1
        ELSEIF(IE.EQ.4) THEN
          IF(SPECTRO.EQ.'PHD') THEN
            IE=1
            I_EXT=2
            IMOD=0
          ELSE
            IE=1
            I_EXT=1
          ENDIF
        ENDIF
      ENDIF
C
      ICALC=IPHI*IE+IPHI*ITHETA+IE*ITHETA
      IF((ICALC.NE.0).AND.(IFTHET.EQ.0)) IRET=3
C
C  When the direction of the analyzer might be experimentally
C    inaccurate, the calculation will be done for nine
C    direction across the one given in the data file
C    with an increment of one degree.
C
      IF(ICHKDIR.EQ.1) THEN
        IF((ITHETA.EQ.1).AND.(IPHI.EQ.0)) THEN
          NPHI=9
          PHI0=PHI0-4.
          PHI1=PHI0+8.
        ELSEIF((IPHI.EQ.1).AND.(ITHETA.EQ.0)) THEN
          NTHETA=9
          THETA0=THETA0-4.
          THETA1=THETA0+8.
        ENDIF
      ENDIF
C
C  Initialization of the values for the scanned angle and the "fixed" one
C
      IF(IPHI.EQ.1) THEN
        N_FIXED=NTHETA
        N_SCAN=NPHI
        FIX0=THETA0
        FIX1=THETA1
        SCAN0=PHI0
        SCAN1=PHI1
        IPH_1=0
      ELSEIF(ITHETA.EQ.1) THEN
        N_FIXED=NPHI
        N_SCAN=NTHETA
        FIX0=PHI0
        FIX1=PHI1
        SCAN0=THETA0
        SCAN1=THETA1
        IPH_1=1
      ELSEIF(IE.EQ.1) THEN
        IF(NTHETA.GE.NPHI) THEN
          N_FIXED=NPHI
          N_SCAN=NTHETA
          FIX0=PHI0
          FIX1=PHI1
          SCAN0=THETA0
          SCAN1=THETA1
          IPH_1=1
        ELSE
          N_FIXED=NTHETA
          N_SCAN=NPHI
          FIX0=THETA0
          FIX1=THETA1
          SCAN0=PHI0
          SCAN1=PHI1
          IPH_1=0
        ENDIF
      ENDIF
C
 602  IF(SPECTRO.EQ.'XAS') THEN
        READ(ICOM,1) RIEN
        READ(ICOM,2) TEXTE5
        READ(ICOM,1) RIEN
C
        READ(ICOM,39) EDGE,NEDGE,INITL,THLUM,PHILUM
        READ(ICOM,19) NE_X,EK_INI,EK_FIN,EPH_INI
C
        LI=NEDGE/2
        IF(NEDGE.GT.1) I_SO=2
        IF(EDGE.EQ.'K') THEN
          NI=1
        ELSEIF(EDGE.EQ.'L') THEN
          NI=2
        ELSEIF(EDGE.EQ.'M') THEN
          NI=3
        ELSEIF(EDGE.EQ.'N') THEN
          NI=4
        ELSEIF(EDGE.EQ.'O') THEN
          NI=5
        ELSEIF(EDGE.EQ.'P') THEN
          NI=6
        ENDIF
      ELSE
        DO JLINE=1,5
          READ(ICOM,1) RIEN
        ENDDO
      ENDIF
C
 603  IF((SPECTRO.EQ.'AED').OR.(SPECTRO.EQ.'APC')) THEN
C
        READ(ICOM,1) RIEN
        READ(ICOM,2) TEXTE6
        READ(ICOM,1) RIEN
C
        READ(ICOM,40) EDGE_C,NEDGE_C,EDGE_I,NEDGE_I,EDGE_A,NEDGE_A
        READ(ICOM,42) I_MULT,IM1,MULT,IM2
        READ(ICOM,14) IPHI_A,ITHETA_A,IFTHET_A,I_INT
        READ(ICOM,15) NPHI_A,NTHETA_A,NFTHET_A
        READ(ICOM,41) PHI0_A,THETA0_A,R1_A
        READ(ICOM,41) PHI1_A,THETA1_A,R2_A
        READ(ICOM,18) IMOD_A,IMOY_A,ACCEPT_A,ICHKDIR_A
C
        LI_C=NEDGE_C/2
        LI_I=NEDGE_I/2
        LI_A=NEDGE_A/2
C
        IF((EDGE_I.EQ.EDGE_A).AND.(LI_I.EQ.LI_A)) THEN
          I_SHELL=1
        ELSE
          I_SHELL=0
        ENDIF
C
        IE_A=0
        NE_A=1
        I_CP_A=0
C
        IF(EDGE_C.EQ.'K') THEN
          AUGER=' '//EDGE_C//EDGE_I//CHAR(48+NEDGE_I)//EDGE_A//CHAR(48+N
     &EDGE_A)
        ELSE
          AUGER=EDGE_C//CHAR(48+NEDGE_C)//EDGE_I//CHAR(48+NEDGE_I)//EDGE
     &_A//CHAR(48+NEDGE_A)
        ENDIF
        AUGER1=AUGER
C
        IF(IPHI_A.EQ.-1) THEN
          IPHI_A=1
          I_EXT_A=0
          ICHKDIR_A=0
          STEREO='YES'
          IF(ABS(PHI1_A-PHI0_A).LT.0.0001) THEN
            PHI0_A=0.
            PHI1_A=360.
            NPHI_A=361
          ENDIF
          IF(ABS(THETA1_A-THETA0_A).LT.0.0001) THEN
            THETA0_A=0.
            THETA1_A=88.
            NTHETA_A=89
          ENDIF
        ELSEIF(IPHI_A.EQ.2) THEN
          IPHI_A=1
          I_EXT_A=1
        ELSEIF(IPHI_A.EQ.3) THEN
          IPHI_A=1
          I_EXT_A=-1
        ELSEIF(ITHETA_A.EQ.2) THEN
          ITHETA_A=1
          I_EXT_A=1
        ELSEIF(ITHETA_A.EQ.3) THEN
          ITHETA_A=1
          I_EXT_A=-1
        ENDIF
C
C  Check for the consistency of the data for the two electrons in
C     APECS, in particular when the sample is rotated (IMOD=1)
C
        IF(SPECTRO.EQ.'APC') THEN
           IF((LI_C.NE.LI).OR.(IMOD_A.NE.IMOD)) THEN
             IRET=11
             GOTO 605
           ENDIF
           DTH=THETA1-THETA0
           DTH_A=THETA1_A-THETA0_A
           DPH=PHI1-PHI0
           DPH_A=PHI1_A-PHI0_A
           IF((IMOD_A.EQ.1).AND.(IPHI_A.NE.IPHI)) IRET=13
           IF((IMOD_A.EQ.1).AND.(ITHETA_A.NE.ITHETA)) IRET=13
           IF((IMOD_A.EQ.1).AND.(NPHI_A.NE.NPHI)) IRET=13
           IF((IMOD_A.EQ.1).AND.(NTHETA_A.NE.NTHETA)) IRET=13
           IF((IMOD_A.EQ.1).AND.(DTH_A.NE.DTH)) IRET=13
           IF((IMOD_A.EQ.1).AND.(DPH_A.NE.DPH)) IRET=13
        ENDIF
C
C  When the direction of the analyzer might be experimentally
C    inaccurate, the calculation will be done for nine
C    direction across the one given in the data file
C    with an increment of one degree.
C
        IF(ICHKDIR_A.EQ.1) THEN
          IF((ITHETA_A.EQ.1).AND.(IPHI_A.EQ.0)) THEN
            NPHI_A=9
            PHI0_A=PHI0_A-4.
            PHI1_A=PHI0_A+8.
          ELSEIF((IPHI_A.EQ.1).AND.(ITHETA_A.EQ.0)) THEN
            NTHETA_A=9
            THETA0_A=THETA0_A-4.
            THETA1_A=THETA0_A+8.
          ENDIF
        ENDIF
C
C  Initialization of the values for the scanned angle and the "fixed" one
C
        IF(IPHI_A.EQ.1) THEN
          N_FIXED_A=NTHETA_A
          N_SCAN_A=NPHI_A
          FIX0_A=THETA0_A
          FIX1_A=THETA1_A
          SCAN0_A=PHI0_A
          SCAN1_A=PHI1_A
          IPH_1_A=0
        ELSEIF(ITHETA_A.EQ.1) THEN
          N_FIXED_A=NPHI_A
          N_SCAN_A=NTHETA_A
          FIX0_A=PHI0_A
          FIX1_A=PHI1_A
          SCAN0_A=THETA0_A
          SCAN1_A=THETA1_A
          IPH_1_A=1
        ENDIF
C
      ELSE
        DO JLINE=1,10
          READ(ICOM,1) RIEN
        ENDDO
      ENDIF
C
      IF(SPECTRO.EQ.'XAS') THEN
        I_CP=1
        NE=NE_X
      ELSE
        I_CP=0
      ENDIF
C
 608  IF(SPECTRO.EQ.'EIG') THEN
C
        READ(ICOM,1) RIEN
        READ(ICOM,2) TEXTE6B
        READ(ICOM,1) RIEN
C
        READ(ICOM,43) NE_EIG,E0_EIG,EFIN_EIG,I_DAMP
C
        NE=NE_EIG
        N_LINE_E=INT((FLOAT(NE_EIG)-0.0001)/4.)+1
        N_LAST=4-(4*N_LINE_E-NE_EIG)
C
        IF(N_LINE_E.GT.1) THEN
          DO JLINE=1,N_LINE_E-1
             J=(JLINE-1)*4
             READ(ICOM,7) I_SPECTRUM(J+1),I_SPECTRUM(J+2),I_SPECTRUM(J+3
     &),I_SPECTRUM(J+4)
          ENDDO
        ENDIF
C
        J=4*(N_LINE_E-1)
C
        READ(ICOM,7) (I_SPECTRUM(J+K), K=1,N_LAST)
C
        READ(ICOM,46) I_PWM,METHOD,RACC,EXPO
        READ(ICOM,47) N_MAX,N_ITER,N_TABLE,SHIFT
        READ(ICOM,48) I_XN,I_VA,I_GN,I_WN
        READ(ICOM,49) LEVIN,ALPHAR,BETAR
C
        ACC=DBLE(RACC)
        IF(ABS(I_PWM).LE.2) THEN
          I_ACC=0
          N_ITER=N_MAX
        ELSEIF(I_PWM.EQ.3) THEN
          I_ACC=1
          N_ITER=N_MAX
        ELSEIF(I_PWM.EQ.-3) THEN
          I_ACC=-1
          N_ITER=N_MAX
        ELSEIF(I_PWM.EQ.4) THEN
          I_ACC=2
        ELSEIF(I_PWM.EQ.-4) THEN
          I_ACC=-2
        ENDIF
        IF(N_MAX.LT.N_ITER) N_ITER=N_MAX
C
        ALPHA=DCMPLX(ALPHAR)
        BETA=DCMPLX(BETAR)
C
C
      ELSE
        DO JLINE=1,9
          READ(ICOM,1) RIEN
        ENDDO
C
      ENDIF
C
 609  READ(ICOM,1) RIEN
      READ(ICOM,2) TEXTE7
      READ(ICOM,1) RIEN
C
      READ(ICOM,21) NO,NDIF,ISPHER,I_GR
      READ(ICOM,50) I_REN,N_REN,REN_R,REN_I
C
      IF(ISPHER.EQ.0) THEN
        IDWSPH=0
        NO=0
      ENDIF
      IF(NO.LT.0) NO=8
      NUMAX(1)=NO/2
C
      READ(ICOM,22) ISFLIP,IR_DIA,ITRTL,I_TEST
C
      IF((SPECTRO.EQ.'AED').OR.(SPECTRO.EQ.'APC')) I_TEST_A=I_TEST
      IF(I_TEST.EQ.1) THEN
        IF(INTERACT.EQ.'DIPOLAR') THEN
          INITL=1
          LI=0
          IPOL=1
        ELSEIF(INTERACT.EQ.'COULOMB') THEN
          LI_C=0
          LI_I=0
        ENDIF
      ENDIF
C
      READ(ICOM,23) NEMET
C
      BACKSPACE ICOM
      NLG=INT((NEMET-0.0001)/3) +1
      DO N=1,NLG
        NRL=3*N
        JD=3*(N-1)+1
        IF(N.EQ.NLG) NRL=NEMET
        READ(ICOM,24) NEMO,(IEMET(J), J=JD, NRL)
        IF(N.EQ.1) NEMET1=NEMO
      ENDDO
C
      READ(ICOM,25) ISOM,NONVOL(JFICH),NPATHP,VINT
C
      IF(I_TEST.EQ.2) THEN
        IF(ABS(IPOL).EQ.1) THEN
          THLUM=-90.
          PHILUM=0.
        ELSEIF(ABS(IPOL).EQ.2) THEN
          THLUM=0.
          PHILUM=0.
        ENDIF
        IMOD=0
        VINT=0.
        A=1.
      ENDIF
C
      IF((NFICHLEC.EQ.1).OR.(IBAS.EQ.1)) ISOM=0
C
      READ(ICOM,26) IFWD,NTHOUT,I_NO,I_RA
C
      IF(NTHOUT.EQ.NDIF-1) IFWD=0
C
      IF(I_RA.EQ.1) NO=0
      DO JAT=1,NAT
        READ(ICOM,27) N_RA(JAT),THFWD(JAT),IBWD(JAT),THBWD(JAT)
        IF(I_RA.EQ.0) THEN
          N_RA(JAT)=NO
          NUMAX(JAT)=NO/2
        ELSEIF(I_RA.EQ.1) THEN
          NUMAX(JAT)=N_RA(JAT)/2
          NO=MAX(N_RA(JAT),NO)
        ENDIF
      ENDDO
C
      READ(ICOM,5) TESLEC
      IF(TESLEC.EQ.'IPW,NCU') THEN
        BACKSPACE ICOM
      ELSE
        IRET=8
        GOTO 605
      ENDIF
C
      READ(ICOM,28) IPW,NCUT,PCTINT,IPP
      READ(ICOM,29) ILENGTH,RLENGTH,UNLENGTH
      READ(ICOM,30) IDWSPH,ISPEED,IATTS,IPRINT
C
      IF(IDWSPH.EQ.0) ISPEED=1
C
      READ(ICOM,31) IDCM,TD,T,RSJ
      READ(ICOM,32) ILPM,XLPM0
C
      IF((IDCM.GE.1).OR.(ILPM.EQ.1)) THEN
        CALL ATDATA
      ENDIF
      NLEC=INT((NAT-0.0001)/4)+1
C
      DO I=1,NLEC
        NDEB=4*(I-1) + 1
        NFIN=MIN0(4*I,NAT)
        READ(ICOM,33) (UJ2(J),J=NDEB,NFIN)
      ENDDO
C
      DO JLINE=1,5
        READ(ICOM,1) RIEN
      ENDDO
      READ(ICOM,5) TESLEC
      IF(TESLEC.EQ.'DATA FI') THEN
        BACKSPACE ICOM
      ELSE
        IRET=9
        GOTO 605
      ENDIF
C
      READ(ICOM,34) INFILE1,IUI1
      READ(ICOM,34) INFILE2,IUI2
      READ(ICOM,34) INFILE3,IUI3
      READ(ICOM,34) INFILE4,IUI4
      READ(ICOM,34) INFILE5,IUI5
      READ(ICOM,34) INFILE6,IUI6
C
      IF(SPECTRO.NE.'APC') THEN
        DO JLINE=1,9
          READ(ICOM,1) RIEN
        ENDDO
      ELSE
        DO JLINE=1,6
          READ(ICOM,1) RIEN
        ENDDO
        READ(ICOM,34) INFILE7,IUI7
        READ(ICOM,34) INFILE8,IUI8
        READ(ICOM,34) INFILE9,IUI9
      ENDIF
C
C  Set up of the switch controlling external
C    reading of the detector directions and
C    averaging over them for an undetected electron
C
      IF(SPECTRO.EQ.'APC') THEN
        IF((I_EXT.EQ.-1).OR.(I_EXT_A.EQ.-1)) THEN
          IF(I_EXT*I_EXT_A.EQ.0) THEN
            WRITE(IUO1,523)
            I_EXT=-1
            I_EXT_A=-1
            OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
            OPEN(UNIT=IUI9, FILE=INFILE9, STATUS='OLD')
            READ(IUI6,713) IDIR,NSET
            READ(IUI9,713) IDIR_A,NSET_A
            IF(IDIR.EQ.2) THEN
              IF(NSET.NE.NSET_A) WRITE(IUO1,524) NSET,NSET_A
              STOP
            ENDIF
          ENDIF
        ENDIF
        IF(I_INT.EQ.1) THEN
          I_EXT=2
        ELSEIF(I_INT.EQ.2) THEN
          I_EXT_A=2
        ELSEIF(I_INT.EQ.3) THEN
          I_EXT=2
          I_EXT_A=2
        ENDIF
      ENDIF
C
      IF(I_EXT.EQ.-1) THEN
        OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
        READ(IUI6,701) IDIR,I_SET,N_POINTS
        READ(IUI6,702) I_PH,N_FIXED,N_SCAN
        DO JS=1,I_SET
          READ(IUI6,703) TH_0(JS),PH_0(JS)
        ENDDO
        CLOSE(IUI6)
        IF(IDIR.NE.2) IRET=12
        IF(I_PH.NE.IPH_1) IPH_1=I_PH
        IF((SPECTRO.EQ.'PHD').OR.(SPECTRO.EQ.'APC')) THEN
          IF(I_PH.EQ.0) THEN
            NTHETA=N_FIXED
            NPHI=N_SCAN
          ELSE
            NTHETA=N_SCAN
            NPHI=N_FIXED
          ENDIF
          ICHKDIR=2
        ENDIF
      ENDIF
      IF(I_EXT.GE.1) THEN
        OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
        READ(IUI6,701) IDIR,I_SET,N_POINTS
        CLOSE(IUI6)
        IF((IDIR.NE.1).AND.(I_EXT.EQ.2)) IRET=12
        N_FIXED=N_POINTS
        N_SCAN=1
        NTHETA=N_POINTS
        NPHI=1
      ENDIF
      IF(I_EXT_A.GE.1) THEN
        IF(SPECTRO.EQ.'APC') THEN
          OPEN(UNIT=IUI9, FILE=INFILE9, STATUS='OLD')
          READ(IUI9,701) IDIR_A,I_SET_A,N_POINTS_A
          CLOSE(IUI9)
        ELSE
          OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
          READ(IUI6,701) IDIR_A,I_SET_A,N_POINTS_A
          CLOSE(IUI6)
        ENDIF
        IF((IDIR_A.NE.1).AND.(I_EXT_A.EQ.2)) IRET=12
        N_FIXED_A=N_POINTS_A
        N_SCAN_A=1
        NTHETA_A=N_POINTS_A
        NPHI_A=1
      ENDIF
C
      IF(I_EXT_A.EQ.-1) THEN
        IF(SPECTRO.EQ.'APC') THEN
          OPEN(UNIT=IUI9, FILE=INFILE9, STATUS='OLD')
          READ(IUI9,701) IDIR_A,I_SET_A,N_POINTS_A
          READ(IUI9,702) I_PH_A,N_FIXED_A,N_SCAN_A
          CLOSE(IUI9)
        ELSE
          OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
          READ(IUI6,701) IDIR_A,I_SET_A,N_POINTS_A
          READ(IUI6,702) I_PH_A,N_FIXED_A,N_SCAN_A
          CLOSE(IUI6)
        ENDIF
        IF(IDIR_A.NE.2) IRET=12
        IF(I_PH_A.EQ.0) THEN
          NTHETA_A=N_FIXED_A
          NPHI_A=N_SCAN_A
        ELSE
          NTHETA_A=N_SCAN_A
          NPHI_A=N_FIXED_A
        ENDIF
        ICHKDIR_A=2
      ENDIF
C
      DO JLINE=1,5
        READ(ICOM,1) RIEN
      ENDDO
C
      READ(ICOM,34) OUTFILE1,IUO1
      READ(ICOM,34) OUTFILE2,IUO2
      READ(ICOM,34) OUTFILE3,IUO3
      READ(ICOM,34) OUTFILE4,IUO4
C
      IUSCR=MAX0(ICOM,IUI2,IUI3,IUI4,IUI5,IUI6,IUI7,IUI8,IUI9,IUO1,IUO2,
     &IUO3,IUO4)+1
      IUSCR2=IUSCR+1
C
      IF(IADS.GE.1) THEN
        OPEN(UNIT=IUI5, FILE=INFILE5, STATUS='OLD')
        READ(IUI5,1) RIEN
        READ(IUI5,12) NATA,NADS1,NADS2,NADS3
        IF(NATA.EQ.1) THEN
          NADS2=0
          NADS3=0
        ELSEIF(NATA.EQ.2) THEN
          NADS3=0
        ENDIF
        READ(IUI5,35) (NZAT(I),I=NAT+1,NAT+NATA)
        READ(IUI5,36) (CHEM(I),I=NAT+1,NAT+NATA)
        READ(IUI5,37) (UJ2(NAT+J),J=1,NATA)
        READ(IUI5,38) NRELA,(PCRELA(I),I=1,NRELA)
        IF(NRELA.EQ.0) THEN
          DO JRELA=1,3
            PCRELA(JRELA)=0.
          ENDDO
        ENDIF
        NADS=NADS1+NADS2+NADS3
        DO JADS=1,NADS
          READ(IUI5,9) (ADS(I,JADS),I=1,3)
        ENDDO
        CLOSE(IUI5)
      ELSE
        NATA=0
        NRELA=0
      ENDIF
C
      GOTO 601
C
 605  REWIND ICOM
      DO JLINE=1,500
        READ(ICOM,5) TESLEC
        IF(TESLEC.EQ.'CONTROL') THEN
          BACKSPACE ICOM
          READ(ICOM,34) OUTFILE1,IUO1
          GOTO 601
        ENDIF
      ENDDO
C
 601  IF((JFICH.EQ.1).OR.(ISOM.EQ.0)) THEN
c     LINE REMOVED BY PYMSSPEC
      ENDIF
      IF((NFICHLEC.GT.1).AND.(ISOM.NE.0)) THEN
         WRITE(IUO1,105) INDATA(JFICH)
      ENDIF
C
      IF(IRET.EQ.1) RETURN 1
      IF(IRET.EQ.3) RETURN 3
      IF(IRET.EQ.4) RETURN 4
      IF(IRET.EQ.5) RETURN 5
      IF(IRET.EQ.6) RETURN 6
      IF(IRET.EQ.7) RETURN 7
      IF(IRET.EQ.8) RETURN 8
      IF(IRET.EQ.9) RETURN 9
      IF(IRET.EQ.10) RETURN 10
      IF(IRET.EQ.11) RETURN 11
      IF(IRET.EQ.12) RETURN 12
      IF(IRET.EQ.13) RETURN 13
C
C
C..........   Writing of the input data in unit IUO1   ..........
C
C
      WRITE(IUO1,100)
      WRITE(IUO1,101)
      WRITE(IUO1,101)
      WRITE(IUO1,102) TEXTE1
      WRITE(IUO1,101)
      WRITE(IUO1,101)
      WRITE(IUO1,203)
C
      IF(I_TEST.NE.2) THEN
        WRITE(IUO1,201) TEXTE2
      ELSE
        IF(ABS(IPOL).EQ.1) THEN
          WRITE(IUO1,525)
        ELSEIF(ABS(IPOL).EQ.2) THEN
          WRITE(IUO1,526)
        ENDIF
      ENDIF
C
      IF(NAT.GT.NATP_M) RETURN 2
      IF(NE.GT.NE_M) RETURN 2
      IF(NEMET.GT.NEMET_M) RETURN 2
C
      IF(I_TEST.EQ.2) GOTO 606
      IF(IBAS.EQ.0) THEN
        WRITE(IUO1,204) A,IBAS
        GOTO 604
      ENDIF
      WRITE(IUO1,103) CRIST,CENTR,IBAS,NAT
      IF(NCRIST.EQ.1) THEN
        BSURA=1.
        CSURA=1.
        WRITE(IUO1,304) A
      ELSEIF((NCRIST.EQ.2).OR.(NCRIST.EQ.7).OR.(NCRIST.EQ.6)) THEN
        BSURA=1.
        WRITE(IUO1,404) A,CSURA
        IF((NCRIST.EQ.6).AND.(CSURA.EQ.1.)) THEN
          WRITE(IUO1,206) ALPHAD
        ELSEIF(NCRIST.EQ.4) THEN
          WRITE(IUO1,306) BETAD
        ENDIF
      ELSEIF((NCRIST.EQ.3).OR.(NCRIST.EQ.5).OR.(NCRIST.EQ.8)) THEN
        WRITE(IUO1,104) A,BSURA,CSURA
        IF(NCRIST.NE.3) THEN
          WRITE(IUO1,106) ALPHAD,BETAD,GAMMAD
        ENDIF
      ENDIF
      IF(NCRIST.EQ.7) THEN
        WRITE(IUO1,107) IH,IK,II,IL
      ELSE
        WRITE(IUO1,207) IH,IK,IL
      ENDIF
      WRITE(IUO1,108) NIV,COUPUR,ITEST,IESURF
      IF(NAT.GT.1) THEN
        DO I=1,NAT
          J=3*(I-1)
          WRITE(IUO1,109) ATBAS(1+J),ATBAS(2+J),ATBAS(3+J),CHEM(I),NZAT(
     &I)
        ENDDO
      ENDIF
      IF(NCRIST.EQ.8) THEN
        DO I=1,3
          J=3*(I-1)
          WRITE(IUO1,209) VECBAS(1+J),VECBAS(2+J),VECBAS(3+J)
        ENDDO
      ENDIF
      IF(IREL.GE.1) THEN
        WRITE(IUO1,110) IREL,NREL,(PCREL(I),I=1,2)
        IF(NREL.GT.2) THEN
          NLIGNE=INT(FLOAT(NREL-2)/4.)+1
          DO J=1,NLIGNE
            WRITE(IUO1,210) (PCREL(I),I=1,4)
          ENDDO
        ENDIF
        IF(NREL.GT.10) RETURN 4
        WRITE(IUO1,112) OMEGAD1,OMEGAD2,IADS
      ENDIF
      IF((IREL.EQ.0).AND.(IADS.EQ.1)) WRITE(IUO1,212) IADS
      IF(IADS.GE.1) THEN
        WRITE(IUO1,501)
        DO JADS=1,NADS
          IF(JADS.LE.NADS1) THEN
            IF(JADS.EQ.1) WRITE(IUO1,303) NAT+1
            WRITE(IUO1,309) (ADS(I,JADS),I=1,3)
          ELSEIF((JADS.GT.NADS1).AND.(JADS.LE.(NADS1+NADS2))) THEN
            IF(JADS.EQ.(NADS1+1)) WRITE(IUO1,303) NAT+2
            WRITE(IUO1,309) (ADS(I,JADS),I=1,3)
          ELSEIF(JADS.GT.(NADS1+NADS2)) THEN
            IF(JADS.EQ.(NADS2+1)) WRITE(IUO1,303) NAT+3
            WRITE(IUO1,309) (ADS(I,JADS),I=1,3)
          ENDIF
        ENDDO
      ENDIF
      IF((IREL.GT.0).OR.(NRELA.GT.0)) WRITE(IUO1,502)
      IF(NRELA.GT.0) THEN
        WRITE(IUO1,311) (PCRELA(I),I=1,NRELA)
      ENDIF
 604  IF(IREL.GT.0) THEN
        WRITE(IUO1,211) (PCREL(I),I=1,NREL)
      ENDIF
C
 606  IF(SPECTRO.EQ.'APC') WRITE(IUO1,517)
C
      IF(SPECTRO.EQ.'PHD') THEN
C
        IF(IPHI.EQ.1) THEN
          IF(STEREO.EQ.' NO') THEN
            WRITE(IUO1,503)
          ELSE
            WRITE(IUO1,527)
          ENDIF
        ENDIF
        IF(IE.EQ.1) WRITE(IUO1,504)
        IF(ITHETA.EQ.1) WRITE(IUO1,505)
        IF(IFTHET.EQ.1) WRITE(IUO1,506)
        IF(I_AMP.EQ.1) WRITE(IUO1,534)
C
        WRITE(IUO1,201) TEXTE4
        WRITE(IUO1,113) ISPIN,IDICHR,IPOL
        WRITE(IUO1,120) NI,NLI,S_O,INITL,I_SO
        WRITE(IUO1,114) IPHI,ITHETA,IE,IFTHET
        WRITE(IUO1,115) NPHI,NTHETA,NE,NFTHET
C
        IF((ITHETA.EQ.1).AND.(IFTHET.EQ.0)) THEN
          IF((THETA0.LT.-90.0).OR.(THETA1.GT.90.0)) THEN
            WRITE(IUO1,508)
            STOP
          ENDIF
          IF(ABS(THLUM).GT.90.0) THEN
            WRITE(IUO1,509)
            STOP
          ENDIF
        ENDIF
C
        WRITE(IUO1,116) PHI0,THETA0,E0,R1
        WRITE(IUO1,216) PHI1,THETA1,EFIN,R2
        WRITE(IUO1,117) THLUM,PHILUM,ELUM
        WRITE(IUO1,118) IMOD,IMOY,ACCEPT,ICHKDIR
C
        IF(IMOY.GT.3) IMOY=3
        IF(IMOY.LT.0) IMOY=0
        IF(IMOY.EQ.0) NDIR=1
        IF(IMOY.EQ.1) NDIR=5
        IF(IMOY.EQ.2) NDIR=13
        IF(IMOY.EQ.3) NDIR=49
        IF((LI.EQ.0).AND.(INITL.NE.0)) INITL=1
C
      ELSEIF(SPECTRO.EQ.'LED') THEN
C
        IF(IPHI.EQ.1) THEN
          IF(STEREO.EQ.' NO') THEN
            WRITE(IUO1,529)
          ELSE
            WRITE(IUO1,530)
          ENDIF
        ENDIF
        IF(IE.EQ.1) WRITE(IUO1,531)
        IF(ITHETA.EQ.1) WRITE(IUO1,532)
        IF(IFTHET.EQ.1) WRITE(IUO1,506)
        IF(I_AMP.EQ.1) WRITE(IUO1,534)
C
        WRITE(IUO1,201) TEXTE4
        WRITE(IUO1,141) ISPIN
        WRITE(IUO1,114) IPHI,ITHETA,IE,IFTHET
        WRITE(IUO1,115) NPHI,NTHETA,NE,NFTHET
C
        IF((ITHETA.EQ.1).AND.(IFTHET.EQ.0)) THEN
          IF((THETA0.LT.-90.0).OR.(THETA1.GT.90.0)) THEN
            WRITE(IUO1,508)
            STOP
          ENDIF
        ENDIF
C
        WRITE(IUO1,116) PHI0,THETA0,E0,R1
        WRITE(IUO1,216) PHI1,THETA1,EFIN,R2
        WRITE(IUO1,142) TH_INI,PHI_INI
        WRITE(IUO1,118) IMOD,IMOY,ACCEPT,ICHKDIR
C
        IF(IMOY.GT.3) IMOY=3
        IF(IMOY.LT.0) IMOY=0
        IF(IMOY.EQ.0) NDIR=1
        IF(IMOY.EQ.1) NDIR=5
        IF(IMOY.EQ.2) NDIR=13
        IF(IMOY.EQ.3) NDIR=49
C
      ELSEIF(SPECTRO.EQ.'XAS') THEN
C
        WRITE(IUO1,507)
        IF(I_AMP.EQ.1) WRITE(IUO1,534)
        WRITE(IUO1,201) TEXTE5
        WRITE(IUO1,113) ISPIN,IDICHR,IPOL
        WRITE(IUO1,134) EDGE,NEDGE,INITL,THLUM,PHILUM
        WRITE(IUO1,119) NE_X,EK_INI,EK_FIN,EPH_INI
C
      ELSEIF(SPECTRO.EQ.'AED') THEN
C
        IF(IPHI_A.EQ.1) THEN
          IF(STEREO.EQ.' NO') THEN
            WRITE(IUO1,515)
          ELSE
            WRITE(IUO1,528)
          ENDIF
        ENDIF
        IF(ITHETA_A.EQ.1) WRITE(IUO1,516)
        IF(I_AMP.EQ.1) WRITE(IUO1,534)
        WRITE(IUO1,201) TEXTE6
        WRITE(IUO1,113) ISPIN,IDICHR,IPOL
        WRITE(IUO1,135) EDGE_C,NEDGE_C,EDGE_I,NEDGE_I,EDGE_A,NEDGE_A
        WRITE(IUO1,140) I_MULT,IM1,MULT,IM2
        WRITE(IUO1,136) IPHI_A,ITHETA_A,IFTHET_A,I_INT
        WRITE(IUO1,137) NPHI_A,NTHETA_A,NFTHET_A
        WRITE(IUO1,138) PHI0_A,THETA0_A,R1_A
        WRITE(IUO1,139) PHI1_A,THETA1_A,R2_A
        WRITE(IUO1,118) IMOD_A,IMOY_A,ACCEPT_A,ICHKDIR_A
C
        IF(IMOY_A.GT.3) IMOY_A=3
        IF(IMOY_A.LT.0) IMOY_A=0
        IF(IMOY_A.EQ.0) NDIR_A=1
        IF(IMOY_A.EQ.1) NDIR_A=5
        IF(IMOY_A.EQ.2) NDIR_A=13
        IF(IMOY_A.EQ.3) NDIR_A=49
C
      ELSEIF(SPECTRO.EQ.'APC') THEN
C
        WRITE(IUO1,518)
        IF(IPHI.EQ.1) WRITE(IUO1,503)
        IF(ITHETA.EQ.1) WRITE(IUO1,505)
        IF(IFTHET.EQ.1) WRITE(IUO1,506)
        IF(I_AMP.EQ.1) WRITE(IUO1,534)
C
        WRITE(IUO1,201) TEXTE4
        WRITE(IUO1,113) ISPIN,IDICHR,IPOL
        WRITE(IUO1,120) NI,NLI,S_O,INITL,I_SO
        WRITE(IUO1,114) IPHI,ITHETA,IE,IFTHET
        WRITE(IUO1,115) NPHI,NTHETA,NE,NFTHET
C
        IF((ITHETA.EQ.1).AND.(IFTHET.EQ.0)) THEN
          IF((THETA0.LT.-90.0).OR.(THETA1.GT.90.0)) THEN
            WRITE(IUO1,508)
            STOP
          ENDIF
          IF(ABS(THLUM).GT.90.0) THEN
            WRITE(IUO1,509)
            STOP
          ENDIF
        ENDIF
C
        WRITE(IUO1,116) PHI0,THETA0,E0,R1
        WRITE(IUO1,216) PHI1,THETA1,EFIN,R2
        WRITE(IUO1,117) THLUM,PHILUM,ELUM
        WRITE(IUO1,118) IMOD,IMOY,ACCEPT,ICHKDIR
C
        IF(IMOY.GT.3) IMOY=3
        IF(IMOY.LT.0) IMOY=0
        IF(IMOY.EQ.0) NDIR=1
        IF(IMOY.EQ.1) NDIR=5
        IF(IMOY.EQ.2) NDIR=13
        IF(IMOY.EQ.3) NDIR=49
        IF((LI.EQ.0).AND.(INITL.NE.0)) INITL=1
C
        WRITE(IUO1,519)
        IF(IPHI_A.EQ.1) WRITE(IUO1,515)
        IF(ITHETA_A.EQ.1) WRITE(IUO1,516)
        WRITE(IUO1,201) TEXTE6
        WRITE(IUO1,113) ISPIN,IDICHR,IPOL
        WRITE(IUO1,135) EDGE_C,NEDGE_C,EDGE_I,NEDGE_I,EDGE_A,NEDGE_A
        WRITE(IUO1,140) I_MULT,IM1,MULT,IM2
        WRITE(IUO1,136) IPHI_A,ITHETA_A,IFTHET_A,I_INT
        WRITE(IUO1,137) NPHI_A,NTHETA_A,NFTHET_A
        WRITE(IUO1,138) PHI0_A,THETA0_A,R1_A
        WRITE(IUO1,139) PHI1_A,THETA1_A,R2_A
        WRITE(IUO1,118) IMOD_A,IMOY_A,ACCEPT_A,ICHKDIR_A
C
        IF(IMOY_A.GT.3) IMOY_A=3
        IF(IMOY_A.LT.0) IMOY_A=0
        IF(IMOY_A.EQ.0) NDIR_A=1
        IF(IMOY_A.EQ.1) NDIR_A=5
        IF(IMOY_A.EQ.2) NDIR_A=13
        IF(IMOY_A.EQ.3) NDIR_A=49
C
        WRITE(IUO1,520)
C
      ELSEIF(SPECTRO.EQ.'EIG') THEN
C
        WRITE(IUO1,143) NE_EIG,E0_EIG,EFIN_EIG,I_DAMP
        DO JLINE=1,N_LINE_E-1
          J=(JLINE-1)*4
          WRITE(IUO1,145) I_SPECTRUM(J+1),I_SPECTRUM(J+2),I_SPECTRUM(J+3
     &),I_SPECTRUM(J+4)
        ENDDO
        J=4*(N_LINE_E-1)
        WRITE(IUO1,145) (I_SPECTRUM(J+K),K=1,N_LAST)
C
        WRITE(IUO1,146) I_PWM,METHOD,RACC,EXPO
        WRITE(IUO1,147) N_MAX,N_ITER,N_TABLE,SHIFT
        WRITE(IUO1,148) I_XN,I_VA,I_GN,I_WN
        WRITE(IUO1,149) LEVIN,ALPHAR,BETAR
        WRITE(IUO1,533)
C
      ENDIF
C
      WRITE(IUO1,201) TEXTE7
C
      IF(SPECTRO.NE.'EIG') THEN
C
        WRITE(IUO1,121) NO,NDIF,ISPHER,I_GR
        WRITE(IUO1,150) I_REN,N_REN,REN_R,REN_I
C
        IF(SPECTRO.EQ.'XAS') NDIF=NDIF+1
C
        WRITE(IUO1,122) ISFLIP,IR_DIA,ITRTL,I_TEST
C
        IF(ISFLIP.EQ.0) THEN
          NSTEP=3
        ELSE
          NSTEP=1
        ENDIF
        DO N=1,NLG
          NRL=3*N
          JD=3*(N-1)+1
          IF(N.EQ.NLG) NRL=NEMET
          IF(N.EQ.1) NEMO=NEMET1
          IF(N.LT.NLG) THEN
          WRITE(IUO1,123) NEMO,(IEMET(J), J=JD, NRL)
          ELSE
            NTE=NEMET-JD+1
            IF(NTE.EQ.1) WRITE(IUO1,223) NEMO,(IEMET(J),J=JD,NEMET)
            IF(NTE.EQ.2) WRITE(IUO1,323) NEMO,(IEMET(J),J=JD,NEMET)
            IF(NTE.EQ.3) WRITE(IUO1,123) NEMO,(IEMET(J),J=JD,NEMET)
          ENDIF
        ENDDO
      ENDIF
      IF(SPECTRO.NE.'EIG') THEN
        WRITE(IUO1,124) ISOM,NONVOL(JFICH),NPATHP,VINT
        WRITE(IUO1,125) IFWD,NTHOUT,I_NO,I_RA
        DO JAT=1,NAT
          WRITE(IUO1,126) N_RA(JAT),THFWD(JAT),IBWD(JAT),THBWD(JAT)
          RTHFWD(JAT)=THFWD(JAT)*PIS180
          RTHBWD(JAT)=THBWD(JAT)*PIS180
        ENDDO
        WRITE(IUO1,127) IPW,NCUT,PCTINT,IPP
        WRITE(IUO1,128) ILENGTH,RLENGTH,UNLENGTH
        WRITE(IUO1,129) IDWSPH,ISPEED,IATTS,IPRINT
      ELSE
        WRITE(IUO1,144) VINT
      ENDIF
      WRITE(IUO1,130) IDCM,TD,T,RSJ
      WRITE(IUO1,131) ILPM,XLPM0
      DO I=1,NLEC
        NDEB=4*(I-1) + 1
        NFIN=4*I
        IF(I.EQ.NLEC) NFIN=NAT
        NUJ=NFIN-NDEB+1
        IF(NUJ.EQ.1) WRITE(IUO1,132) (UJ2(J),J=NDEB,NFIN)
        IF(NUJ.EQ.2) WRITE(IUO1,232) (UJ2(J),J=NDEB,NFIN)
        IF(NUJ.EQ.3) WRITE(IUO1,332) (UJ2(J),J=NDEB,NFIN)
        IF(NUJ.EQ.4) WRITE(IUO1,432) (UJ2(J),J=NDEB,NFIN)
      ENDDO
      IF(IADS.EQ.1) THEN
        IF(NATA.EQ.1) WRITE(IUO1,133) (UJ2(J),J=NAT+1,NAT+NATA)
        IF(NATA.EQ.2) WRITE(IUO1,233) (UJ2(J),J=NAT+1,NAT+NATA)
        IF(NATA.EQ.3) WRITE(IUO1,333) (UJ2(J),J=NAT+1,NAT+NATA)
      ENDIF
C
      IF(UNLENGTH.EQ.'ATU') RLENGTH=RLENGTH*BOHR/A
      IF(UNLENGTH.EQ.'ANG') RLENGTH=RLENGTH/A
      IF(IBAS.GT.0) THEN
        OMEGA1=OMEGAD1*PIS180
        OMEGA2=OMEGAD2*PIS180
      ENDIF
      QD=0.
      DO J=1,NATM
        UJ2(J)=UJ2(J)/(A*A)
      ENDDO
      IF(E0.EQ.0.) E0=0.0001
      NPOINT=NPHI*NE*NTHETA
      ISORT1=0
      IF(NPOINT.GT.250) THEN
        ISORT1=1
        WRITE(IUO1,510)
      ENDIF
C
      IF(IDWSPH.EQ.1) THEN
        NFAC=N_GAUNT
      ELSE
        NFAC=4*NL_M
      ENDIF
      IF(SPECTRO.EQ.'EIG') THEN
C
C  Switch for including vibrational damping into the MS matrix
C
C           I_VIB = 0 : no vibrations included
C           I_VIB = 1 : vibrations included
C
C         and mean free path-like damping
C
C           I_MFP = 0 : no Im(k) damping included
C           I_MFP = 1 : Im(k) damping included
C
        I_VIB=MOD(I_DAMP,2)
        IF(I_VIB.EQ.1) THEN
          IDWSPH=1
        ELSE
          IDWSPH=0
        ENDIF
        IF(I_DAMP.LE.1) THEN
          I_MFP=0
        ELSE
          I_MFP=1
        ENDIF
      ENDIF
C
C  Computing the renormalization coefficients
C
      IF(I_REN.LE.4) THEN
        CALL COEF_RENORM(NDIF)
      ELSEIF(I_REN.EQ.5) THEN
        CALL COEF_LOEWDIN(NDIF)
      ELSEIF(I_REN.GE.6.AND.I_REN.LE.8) THEN
        CALL COEF_SIGMAPI2(NDIF)
      ENDIF
C
C  Storage of the logarithm of the Gamma function GLD(N+1,N_INT)
C  for integer (N_INT=1) and semi-integer (N_INT=2) values :
C
C    GLD(N+1,1)   =   Log(N!) for N integer
C    GLD(N+1/2,2) =   Log(N!) for N semi-integer
C
      IF((ISPHER.GE.0).OR.(I_MULT.EQ.1)) THEN
        GLG(1)=0.0
        GLD(1,1)=0.D0
        GLD(1,2)=DLOG(SQPI/2.D0)
        DO I=2,NFAC
          J=I-1
          GLG(I)=GLG(J)+ALOG(FLOAT(J))
          GLD(I,1)=GLD(J,1)+DLOG(DFLOAT(J))
          GLD(I,2)=GLD(J,2)+DLOG(DFLOAT(J) +0.5D0)
        ENDDO
      ELSEIF((IFTHET.EQ.1).AND.(ITEST.EQ.1)) THEN
        GLG(1)=0.0
        DO I=2,NFAC
          J=I-1
          GLG(I)=GLG(J)+ALOG(FLOAT(J))
        ENDDO
      ENDIF
      EXPF(0,0)=1.
      EXPR(0,0)=1.
      FACT1L=0.D0
      DO L=1,2*NL_M-2
        XDEN=1./SQRT(FLOAT(L+L+1))
        DXDEN=1.D0/DSQRT(DFLOAT(L+L+1))
        FACT1L=FACT1L+DLOG(DFLOAT(L))
        FACT2L=DLOG(DFLOAT(L+1))
        DO M1=0,L
          EXPF(M1,L)=EXP(0.5*(GLG(L+M1+1)-GLG(L-M1+1)))
          DEXPF=DEXP(0.5D0*(GLD(L+M1+1,1)-GLD(L-M1+1,1)))
          EXPR(M1,L)=EXP(0.5*(GLG(L+L+1)-GLG(L+M1+1)-GLG(L-M1+1)))
          EXPF2(L,M1)=EXPF(M1,L)*XDEN
          DEXPF2(L,M1)=DEXPF*DXDEN
          IF(M1.GT.0) THEN
            FACT2L=FACT2L+DLOG(DFLOAT(1+L+M1))
          ENDIF
          IF(L.LT.NL_M) THEN
            DO M2=0,L
              CF(L,M1,M2)=SQRT(FLOAT((L*L-M1*M1)*(L*L-M2*M2)))/FLOAT(L)
            ENDDO
          ENDIF
        ENDDO
        FSQ(L)=EXP(0.5*REAL(FACT2L-FACT1L))
        DFSQ(L)=DEXP(0.5D0*(FACT2L-FACT1L))
      ENDDO
C
      IF((INITL.LT.-1).OR.(INITL.GT.2)) THEN
        INITL=1
        WRITE(IUO1,511)
      ENDIF
      NEPS=2-ABS(IPOL)
      IF(IDICHR.GE.1) NEPS=1
      ISTEP_LF=ABS(INITL)
      IF(INITL.EQ.-1) THEN
        LF1=LI-1
        LF2=LF1
      ELSEIF(INITL.EQ.1) THEN
        LF1=LI+1
        LF2=LF1
      ELSEIF(INITL.EQ.2) THEN
        LF1=LI-1
        LF2=LI+1
      ELSEIF(INITL.EQ.0) THEN
        LF1=LI
        LF2=LI
        ISTEP_LF=1
      ENDIF
C
C  Initialization of the values of ji if spin-orbit is taken
C                  into account.
C
C  Here :   JI is the loop index going from JF1 to JF2 with :
C
C                   JI=1    : ji = li + 1/2
C                   JI=2    : ji = li - 1/2
C
      IF(I_SO.EQ.0) THEN
        JF1=1
        JF2=2
      ELSEIF(I_SO.EQ.1) THEN
        IF(S_O.EQ.'1/2') THEN
          IF(LI.EQ.0) THEN
            JF1=1
            JF2=1
          ELSEIF(LI.EQ.1) THEN
            JF1=2
            JF2=2
          ENDIF
        ELSEIF(S_O.EQ.'3/2') THEN
          IF(LI.EQ.1) THEN
            JF1=1
            JF2=1
          ELSEIF(LI.EQ.2) THEN
            JF1=2
            JF2=2
          ENDIF
        ELSEIF(S_O.EQ.'5/2') THEN
          IF(LI.EQ.2) THEN
            JF1=1
            JF2=1
          ELSEIF(LI.EQ.3) THEN
            JF1=2
            JF2=2
          ENDIF
        ELSEIF(S_O.EQ.'7/2') THEN
          IF(LI.EQ.3) THEN
            JF1=1
            JF2=1
          ELSEIF(LI.EQ.4) THEN
            JF1=2
            JF2=2
          ENDIF
        ELSEIF(S_O.EQ.'9/2') THEN
          IF(LI.EQ.4) THEN
            JF1=1
            JF2=1
          ELSE
            RETURN 7
          ENDIF
        ELSE
          RETURN 7
        ENDIF
      ELSEIF(I_SO.EQ.2) THEN
        JF1=1
        JF2=2
      ELSE
        RETURN 7
      ENDIF
C
      IF(NI.LE.5) THEN
         NNL=NI*(NI-1)/2 +LI+1
      ELSEIF(NI.EQ.6) THEN
         NNL=NI*(NI-1)/2 +LI
      ELSEIF(NI.EQ.7) THEN
         NNL=NI*(NI-1)/2 +LI-3
      ENDIF
C
C  Storage of the Clebsch-Gordan coefficients for the spin-orbit
C  dependent coupling matrix elements in the array CG(MJI,JI,JSPIN).
C
C  Here :           JI=1    : ji = li + 1/2
C                   JI=2    : ji = li - 1/2
C                   MJI     : mji + 1/2
C                   JSPIN=1 : msi = +1/2
C                   JSPIN=2 : msi = -1/2
C
C              so that all indices remain integer
C
      IF((I_SO.GT.0).OR.(ISPIN.EQ.1).OR.(SPECTRO.EQ.'APC')) THEN
        DO JS=1,2
          DO JI=1,2
            DO MJI=-LI,LI+1
              CG(MJI,JI,JS)=0.0
            ENDDO
          ENDDO
        ENDDO
        DO MJI=-LI,LI+1
          CG(MJI,1,1)=SQRT(FLOAT(LI+MJI)/FLOAT(LI+LI+1))
          CG(MJI,1,2)=SQRT(FLOAT(LI-MJI+1)/FLOAT(LI+LI+1))
          IF((MJI.GT.-LI).AND.(MJI.LT.LI+1)) THEN
            CG(MJI,2,1)=-SQRT(FLOAT(LI-MJI+1)/FLOAT(LI+LI+1))
            CG(MJI,2,2)=SQRT(FLOAT(LI+MJI)/FLOAT(LI+LI+1))
          ENDIF
        ENDDO
      ENDIF
C
C
C  Storage of the Clebsch-Gordan coefficients for the Auger multiplet
C  dependent coupling matrix elements in the array CGA(LJ1,MJ1,LJ2,MJ2,LJ).
C
C  Here :           LJ1  is an integer index related to J1 (LJ1=2*J1)
C                   LMJ1  is an integer index related to MJ1 (LMJ1=2*MJ1)
C                   LJ2  is an integer index related to J2 (LJ2=2*J2)
C                   LMJ2  is an integer index related to MJ2 (LMJ2=2*MJ2)
C                   LJ is an integer index related to J :
C                         J = FLOAT(LJ) for J integer
C                         J = FLOAT(LJ) + 0.5 for J half integer
C
C              so that all indices remain integer
C
      IF((SPECTRO.EQ.'AED').OR.(SPECTRO.EQ.'APC')) THEN
        IF(I_MULT.EQ.1) THEN
          N=3
          MJ3=0.D0
          LJ_MAX=2*(LI_I+LI_A+1)
          DO LJ1=0,LJ_MAX
            J1=DFLOAT(LJ1)/2.D0
            DO LMJ1=-LJ1,LJ1,2
              MJ1=DFLOAT(LMJ1)/2.D0
              DO LJ2=0,LJ_MAX
                J2=DFLOAT(LJ2)/2.D0
                DO LMJ2=-LJ2,LJ2,2
                  MJ2=DFLOAT(LMJ2)/2.D0
                  CALL N_J(J1,MJ1,J2,MJ2,MJ3,NJ,I_INT,N)
C
                  JJ12=J1-J2
                  JL12=MJ1-MJ2
C
                  LJ12=INT(JJ12+SIGN(SMALL,JJ12))
                  LL12=INT(JL12+SIGN(SMALL,JL12))
C
                  JJ_MIN=ABS(LJ12)
                  JJ_MAX=J1+J2
                  LJJ_MIN=INT(JJ_MIN+SIGN(SMALL,JJ_MIN))
                  LJJ_MAX=INT(JJ_MAX+SIGN(SMALL,JJ_MAX))
C
                  DO LJJ=LJJ_MIN,LJJ_MAX,1
                    IF(I_INT.EQ.1) THEN
                      JJ=DFLOAT(LJJ)
                    ELSE
                      JJ=DFLOAT(LJJ)+0.5D0
                    ENDIF
                    L_EXP=INT(J1-J2+MJ1+MJ2)
                    IF(MOD(L_EXP,2).EQ.0) THEN
                      CGA(LJ1,LMJ1,LJ2,LMJ2,LJJ)=NJ(LJJ)*SQRT(2.*REAL(JJ
     &)+1.)
                    ELSE
                      CGA(LJ1,LMJ1,LJ2,LMJ2,LJJ)=-NJ(LJJ)*SQRT(2.*REAL(J
     &J)+1.)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C  Storage of another of the spin Clebsch-Gordan used
C    when the Auger line is multiplet-resolved. It
C    originates from the coupling of SA and SC,
C    the spins of the Auger electron of the original
C    core electron (which is supposed to be the same
C    as that of the photoelectron).
C
C    CG_S(I,J,K) with : I = 1 ---> MSA = -1/2
C                       I = 2 ---> MSA =  1/2
C                       J = 1 ---> MSC = -1/2
C                       J = 2 ---> MSC =  1/2
C                       K = 1 ---> S   =  0
C                       K = 2 ---> S   =  1
C
C                       MS = MSA+MSC
C
      IF(I_MULT.EQ.1) THEN
        CG_S(1,1,1)=0.
        CG_S(1,1,2)=1.
        CG_S(1,2,1)=-0.707107
        CG_S(1,2,2)= 0.707107
        CG_S(2,1,1)= 0.707107
        CG_S(2,1,2)= 0.707107
        CG_S(2,2,1)= 0.
        CG_S(2,2,2)= 1.
      ENDIF
C
C  Initialization of the variables used when only one multiplet
C    is taken into account in the Auger peak
C
      IF(I_MULT.EQ.1) THEN
        MULTIPLET=CHAR(48+IM1)//MULT//CHAR(48+IM2)
        IF(MOD(IM1,2).EQ.0) THEN
          WRITE(IUO1,522) IM1
          STOP
        ENDIF
        S_MUL=(IM1-1)/2
        J_MUL=IM2
        IF(MULT.EQ.'S') THEN
          L_MUL=0
        ELSEIF(MULT.EQ.'P') THEN
          L_MUL=1
        ELSEIF(MULT.EQ.'D') THEN
          L_MUL=2
        ELSEIF(MULT.EQ.'F') THEN
          L_MUL=3
        ELSEIF(MULT.EQ.'G') THEN
          L_MUL=4
        ELSEIF(MULT.EQ.'H') THEN
          L_MUL=5
        ELSEIF(MULT.EQ.'I') THEN
          L_MUL=6
        ELSEIF(MULT.EQ.'K') THEN
          L_MUL=7
        ELSEIF(MULT.EQ.'L') THEN
          L_MUL=8
        ELSEIF(MULT.EQ.'M') THEN
          L_MUL=9
        ELSE
          WRITE(IUO1,521) MULTIPLET
          STOP
        ENDIF
      ENDIF
C
C..........  Check of the dimensioning in the Gaussian case  ..........
C
      CALL STOP_EXT(I_EXT,I_EXT_A,SPECTRO)
C
C....................   Read FORMAT   ....................
C
C
   1  FORMAT(A7)
   2  FORMAT(21X,10A4)
   3  FORMAT(7X,A3,9X,A1,9X,I1,6X,I4)
   4  FORMAT(8X,F6.3,4X,F6.3,4X,F6.3,3X,A3)
   5  FORMAT(49X,A7)
   6  FORMAT(7X,F6.2,4X,F6.2,4X,F6.2)
   7  FORMAT(8X,I2,8X,I2,8X,I2,8X,I2)
   8  FORMAT(8X,I2,8X,F6.3,3X,I3,9X,I1)
   9  FORMAT(8X,F9.6,1X,F9.6,1X,F9.6,2X,A2,2X,I2)
  10  FORMAT(9X,I1,8X,I2,7X,F5.1,5X,F5.1)
  11  FORMAT(7X,F5.1,3(5X,F5.1))
  12  FORMAT(7X,F6.2,4X,F6.2,6X,I1)
  13  FORMAT(7X,A3,9X,I1,9X,I1,8X,I2)
  14  FORMAT(8X,I2,9X,I1,9X,I1,9X,I1,9X,I1)
  15  FORMAT(7X,I3,7X,I3,7X,I3,7X,I3)
  16  FORMAT(6X,F7.2,3X,F7.2,3X,F7.2,5X,F6.3)
  17  FORMAT(6X,F7.2,3X,F7.2,2X,F8.2)
  18  FORMAT(9X,I1,9X,I1,8X,F5.2,6X,I1)
  19  FORMAT(7X,I3,6X,F7.2,3X,F7.2,2X,F8.2)
  20  FORMAT(8X,I1,A1,8X,A3,7X,I2,8X,I2)
  21  FORMAT(8X,I2,8X,I2,9X,I1,9X,I1)
  22  FORMAT(9X,I1,9X,I1,9X,I1,9X,I1)
  23  FORMAT(8X,I2)
  24  FORMAT(8X,I2,3(8X,I2))
  25  FORMAT(9X,I1,8X,I2,1X,I9,8X,F6.2)
  26  FORMAT(9X,I1,9X,I1,9X,I1,9X,I1)
  27  FORMAT(9X,I1,6X,F6.2,7X,I1,7X,F6.2)
  28  FORMAT(9X,I1,9X,I1,7X,F8.4,4X,I1)
  29  FORMAT(9X,I1,7X,F6.2,4X,A3)
  30  FORMAT(9X,I1,8X,I2,9X,I1,9X,I1)
  31  FORMAT(9X,I1,6X,F8.3,2X,F8.3,5X,F4.2)
  32  FORMAT(8X,I2,7X,F6.2)
  33  FORMAT(8X,F8.5,2X,F8.5,2X,F8.5,2X,F8.5)
  34  FORMAT(9X,A24,5X,I2)
  35  FORMAT(18X,I2,8X,I2,8X,I2)
  36  FORMAT(18X,A2,8X,A2,8X,A2)
  37  FORMAT(18X,F8.5,2X,F8.5,2X,F8.5)
  38  FORMAT(9X,I1,7X,F5.1,5X,F5.1,5X,F5.1)
  39  FORMAT(8X,A1,I1,8X,I2,6X,F7.2,3X,F7.2)
  40  FORMAT(8X,A1,I1,8X,A1,I1,8X,A1,I1)
  41  FORMAT(6X,F7.2,3X,F7.2,5X,F6.3)
  42  FORMAT(9X,I1,8X,I1,A1,I1)
  43  FORMAT(7X,I3,6X,F7.2,3X,F7.2,6X,I1)
  44  FORMAT(9X,I1)
  46  FORMAT(8X,I2,6X,A4,9X,F7.5,2X,F6.3)
  47  FORMAT(5X,I5,6X,I4,6X,I4,8X,F6.3)
  48  FORMAT(9X,I1,9X,I1,9X,I1,9X,I1)
  49  FORMAT(8X,I2,6X,F7.2,3X,F7.2)
  50  FORMAT(9X,I1,9X,I1,6X,F8.3,2X,F8.3)
C
C
C....................   Write FORMAT   ....................
C
C
 100  FORMAT(//////////,'******************************', '*************
     &***************************************')
 101  FORMAT('*********************',40X,'*********************')
 102  FORMAT('*********************',10A4,'*********************')
 103  FORMAT(10X,A3,9X,A1,9X,I1,6X,I4,9X,'CRIST,CENTR,IBAS,NAT')
 104  FORMAT(11X,F6.3,4X,F6.3,4X,F6.3,15X,'A,BSURA,CSURA')
 105  FORMAT(///,'ooooooooooooooooooooooooooooooooooooooooo','oooooooooo
     &ooooooooooooooooooooooooooooooo',/,'oooooooooooooooo',50X,'ooooooo
     &ooooooooo',/,'oooooooooooooooo   INPUT DATA FILE   :  ',A24,'  ooo
     &ooooooooooooo',/,'oooooooooooooooo',50X,'oooooooooooooooo',/,'oooo
     &oooooooooooooooooooooooo','ooooooooooooooooooooooooooooooooooooooo
     &oooooooooo','ooooo',///)
 106  FORMAT(10X,F6.2,4X,F6.2,4X,F6.2,16X,'ALPHAD,BETAD,GAMMAD')
 107  FORMAT(11X,I2,8X,I2,8X,I2,8X,I2,9X,'H,K,I,L')
 108  FORMAT(12X,I1,8X,F6.3,3X,I3,9X,I1,9X,'NIV,COUPUR,ITEST,IESURF')
 109  FORMAT(11X,F9.6,1X,F9.6,1X,F9.6,2X,A2,2X,I2,4X,'ATBAS,CHEM(NAT)','
     &,NZAT(NAT)')
 110  FORMAT(12X,I1,8X,I2,7X,F5.1,5X,F5.1,7X,'IREL,NREL,PCREL(NREL)')
 112  FORMAT(10X,F6.2,4X,F6.2,6X,I1,19X,'OMEGA1,OMEGA2,IADS')
 113  FORMAT(12X,I1,9X,I1,8X,I2,19X,'ISPIN,IDICHR,IPOL')
 114  FORMAT(11X,I2,9X,I1,9X,I1,9X,I1,9X,'IPHI,ITHETA,IE,',             
     &'IFTHET')
 115  FORMAT(10X,I3,7X,I3,7X,I3,7X,I3,9X,'NPHI,NTHETA,NE,NFTHET')
 116  FORMAT(9X,F7.2,3X,F7.2,3X,F7.2,5X,F6.3,5X,'PHI0,THETA0,E0,R0')
 117  FORMAT(9X,F7.2,3X,F7.2,2X,F8.2,16X,'THLUM,PHILUM,ELUM')
 118  FORMAT(12X,I1,9X,I1,8X,F5.2,6X,I1,9X,'IMOD,IMOY,ACCEPT,ICHKDIR')
 119  FORMAT(10X,I3,6X,F7.2,3X,F7.2,2X,F8.2,6X,'NE,EK_INI,','EK_FIN,EPH_
     &INI')
 120  FORMAT(11X,I1,A1,8X,A3,7X,I2,8X,I2,9X,'LI,S-O,INITL,I_SO')
 121  FORMAT(11X,I2,8X,I2,9X,I1,9X,I1,9X,'NO,NDIF,ISPHER,I_GR')
 122  FORMAT(12X,I1,9X,I1,9X,I1,9X,I1,9X,'ISFLIP,IR_DIA,ITRTL,I_TEST')
 123  FORMAT(11X,I2,3(8X,I2),9X,'NEMET,IEMET(NEMET)')
 124  FORMAT(12X,I1,8X,I2,6X,I4,7X,F6.2,6X,'ISOM,NONVOL,NPATH,VINT')
 125  FORMAT(12X,I1,9X,I1,9X,I1,9X,I1,9X,'IFWD,NTHOUT,I_NO,I_RA')
 126  FORMAT(12X,I1,7X,F6.2,6X,I1,7X,F6.2,6X,'N_RA(NAT),THFWD(NAT)',',IB
     &WD(NAT),THBWD(NAT)')
 127  FORMAT(12X,I1,9X,I1,7X,F8.4,4X,I1,9X,'IPW,NCUT,PCTINT,IPP')
 128  FORMAT(12X,I1,7X,F6.2,4X,A3,19X,'ILENGTH,RLENGTH,UNLENGTH')
 129  FORMAT(12X,I1,8X,I2,9X,I1,9X,I1,9X,'IDWSPH,ISPEED,IATT,IPRINT')
 130  FORMAT(12X,I1,6X,F8.3,2X,F8.3,5X,F4.2,6X,'IDCM,TD,T,RSJ')
 131  FORMAT(11X,I2,7X,F6.2,26X,'ILPM,XLPM0')
 132  FORMAT(11X,F8.5,33X,'UJ2(NAT)  : ','SUBSTRATE')
 133  FORMAT(11X,F8.5,33X,'UJ2(NATA) : ','ADSORBATES')
 134  FORMAT(11X,A1,I1,8X,I2,6X,F7.2,3X,F7.2,6X,'EDGE,INITL,THLUM,','PHI
     &LUM')
 135  FORMAT(11X,A1,I1,8X,A1,I1,8X,A1,I1,19X,'EDGE_C,EDGE_I,','EDGE_A')
 136  FORMAT(11X,I2,9X,I1,9X,I1,9X,I1,9X,'IPHI_A,ITHETA_A,','IFTHET_A,I_
     &INT')
 137  FORMAT(10X,I3,7X,I3,7X,I3,19X,'NPHI_A,NTHETA_A,NFTHET_A')
 138  FORMAT(9X,F7.2,3X,F7.2,5X,F6.3,15X,'PHI0_A,THETA0_A,R0_A')
 139  FORMAT(9X,F7.2,3X,F7.2,5X,F6.3,15X,'PHI1_A,THETA1_A,R1_A')
 140  FORMAT(12X,I1,8X,I1,A1,I1,28X,'I_MULT,MULT')
 141  FORMAT(12X,I1,39X,'ISPIN')
 142  FORMAT(9X,F7.2,3X,F7.2,26X,'TH_INI,PHI_INI')
 143  FORMAT(10X,I3,6X,F7.2,3X,F7.2,6X,I1,9X,'NE,EK_INI,EK_FIN,I_DAMP')
 144  FORMAT(10X,F6.2,36X,'VINT')
 145  FORMAT(11X,I2,8X,I2,8X,I2,8X,I2,9X,'I_SPECTRUM(NE)')
 146  FORMAT(11X,I2,6X,A4,9X,F7.5,2X,F6.3,5X,'I_PWM,METHOD,ACC,EXPO')
 147  FORMAT(8X,I5,6X,I4,6X,I4,8X,F6.3,5X,'N_MAX,N_ITER,N_TABLE,SHIFT')
 148  FORMAT(12X,I1,9X,I1,9X,I1,9X,I1,9X,'I_XN,I_VA,I_GN,I_WN')
 149  FORMAT(11X,I2,6X,F7.2,3X,F7.2,16X,'L,ALPHA,BETA')
 150  FORMAT(12X,I1,9X,I1,6X,F8.3,2X,F8.3,5X,'I_REN,N_REN,REN_R,REN_I')
C
 201  FORMAT(///,21X,10A4,////)
 203  FORMAT('**************************************************',      
     &'********************************',//////////)
 204  FORMAT(11X,F6.3,5X,I1,29X,'A,IBAS')
 206  FORMAT(10X,F6.2,36X,'ALPHAD')
 207  FORMAT(11X,I2,8X,I2,8X,I2,19X,'H,K,L')
 209  FORMAT(11X,F9.6,1X,F9.6,1X,F9.6,12X,'VECBAS')
 210  FORMAT(10X,F5.1,3(5X,F5.1),7X,'PCREL(NREL)')
 211  FORMAT(20X,'SUBSTRATE : ',10(F5.1,','))
 212  FORMAT(32X,I1,19X,'IADS')
 216  FORMAT(9X,F7.2,3X,F7.2,3X,F7.2,5X,F6.3,5X,'PHI1,THETA1,EFIN,R1')
 223  FORMAT(11X,I2,1(8X,I2),29X,'NEMET,IEMET(NEMET)')
 232  FORMAT(11X,F8.5,2X,F8.5,23X,'UJ2(NAT)  : ','SUBSTRATE')
 233  FORMAT(11X,F8.5,2X,F8.5,23X,'UJ2(NATA) : ','ADSORBATES')
C
 303  FORMAT(/,33X,'ATOMS OF TYPE ',I1,' :',/)
 304  FORMAT(11X,F6.3,35X,'A')
 306  FORMAT(10X,F6.2,36X,'BETAD')
 309  FORMAT(11X,F9.6,1X,F9.6,1X,F9.6,12X,'XADS,YADS,ZADS')
 311  FORMAT(20X,'ADSORBATE : ',3(F5.1,','))
 323  FORMAT(11X,I2,2(8X,I2),19X,'NEMET,IEMET(NEMET)')
 332  FORMAT(11X,F8.5,2X,F8.5,2X,F8.5,13X,'UJ2(NAT)  : ','SUBSTRATE')
 333  FORMAT(11X,F8.5,2X,F8.5,2X,F8.5,13X,'UJ2(NATA) : ','ADSORBATES')
C
 404  FORMAT(11X,F6.3,4X,F6.3,25X,'A,CSURA')
 432  FORMAT(11X,F8.5,2X,F8.5,2X,F8.5,2X,F8.5,3X,'UJ2(NAT)  : ','SUBSTRA
     &TE')
C
 501  FORMAT(//,30X,'POSITION OF THE ADSORBATES :')
 502  FORMAT(///,25X,'VALUE OF THE RELAXATIONS :',/)
 503  FORMAT(///,14X,'TYPE OF CALCULATION : AZIMUTHAL PHOTOELECTRON',' D
     &IFFRACTION')
 504  FORMAT(///,18X,'TYPE OF CALCULATION : FINE STRUCTURE ','OSCILLATIO
     &NS')
 505  FORMAT(///,16X,'TYPE OF CALCULATION : POLAR PHOTOELECTRON',' DIFFR
     &ACTION')
 506  FORMAT(///,23X,'TYPE OF CALCULATION : SCATTERING FACTOR')
 507  FORMAT(///,28X,'TYPE OF CALCULATION : EXAFS')
 508  FORMAT(///,2X,' <<<<<<<<<<  THE THETA VARIATION EXCEEDS THE ',  'P
     &HYSICAL LIMITS (-90,+90)  >>>>>>>>>>',///)
 509  FORMAT(///,2X,' <<<<<<<<<<  THE THLUM VARIATION EXCEEDS THE ',  'P
     &HYSICAL LIMITS (-90,+90)  >>>>>>>>>>',///)
 510  FORMAT(///,4X,' <<<<<<<<<<  AS THE CALCULATION HAS MORE THAN ','25
     &0 POINTS, SOME OUTPUTS HAVE BEEN SUPRESSED  >>>>>>>>>>',///)
 511  FORMAT(///,4X,' <<<<<<<<<<  INCORRECT VALUE OF INITL, THE ',    'C
     &ALCULATION IS PERFORMED WITH INITL = 1  >>>>>>>>>>')
 512  FORMAT(///,4X,' <<<<<<<<<<  IMPOSSIBLE TO HAVE A SPIN RESOLVED ','
     &EXAFS EXPERIMENT : DECREASE IDICHR  >>>>>>>>>>')
 513  FORMAT(///,15X,' <<<<<<<<<<  IMPOSSIBLE TO HAVE IPOL = 0 AND ','ID
     &ICHR > 0  >>>>>>>>>>')
 514  FORMAT(///,15X,' <<<<<<<<<<  IMPOSSIBLE TO HAVE IDICHR = 2 AND ','
     &ISPIN = 0  >>>>>>>>>>')
 515  FORMAT(///,12X,'TYPE OF CALCULATION : AZIMUTHAL AUGER ELECTRON',' 
     &DIFFRACTION')
 516  FORMAT(///,16X,'TYPE OF CALCULATION : POLAR AUGER ELECTRON',' DIFF
     &RACTION')
 517  FORMAT(///,10X,'TYPE OF CALCULATION : AUGER PHOTOELECTRON ','COINC
     &IDENCE SPECTROSCOPY')
 518  FORMAT(///,9X,'------------------------  FIRST ELECTRON : ','-----
     &-------------------')
 519  FORMAT(///,9X,'------------------------ SECOND ELECTRON : ','-----
     &-------------------')
 520  FORMAT(///,9X,'----------------------------------------------','--
     &--------------------')
 521  FORMAT(///,4X,' <<<<<<<<<<  ',A3,' IS NOT IMPLEMENTED IN THIS ','V
     &ERSION  >>>>>>>>>>')
 522  FORMAT(///,4X,' <<<<<<<<<<  WRONG NAME FOR THE MULTIPLET','  >>>>>
     &>>>>>',/,4X,' <<<<<<<<<<  ODD NUMBER ','EXPECTED INSTEAD OF',I2,' 
     & >>>>>>>>>>')
 523  FORMAT(///,4X,' <<<<<<<<<<  BOTH DETECTOR DIRECTIONS MUST BE ','EI
     &THER INTERNAL OR EXTERNAL  >>>>>>>>>>',/,8X,' -----> PROCEEDING WI
     &TH EXTERNAL DIRECTIONS',/)
 524  FORMAT(///,4X,' <<<<<<<<<<  AVERAGING OVER ',I3,' DOMAINS ','FOR P
     &HOTOELECTRON   >>>>>>>>>>',/,4X,' <<<<<<<<<<  AVERAGING OVER ',I3,
     &' DOMAINS ','FOR AUGER ELECTRON   >>>>>>>>>>',/,8X,' -----> IMPOSS
     &IBLE : CHECK INPUT FILES !')
 525  FORMAT(///,14X,'ATOMIC CALCULATION : Z AXIS ALONG POLARIZATION ','
     &DIRECTION',/,'  ',/,'   ',/,'  ')
 526  FORMAT(///,18X,'ATOMIC CALCULATION : Z AXIS ALONG LIGHT ','DIRECTI
     &ON',/,'  ',/,'   ',/,'  ')
 527  FORMAT(///,11X,'TYPE OF CALCULATION : FULL HEMISPHERE',' PHOTOELEC
     &TRON DIFFRACTION')
 528  FORMAT(///,10X,'TYPE OF CALCULATION : FULL HEMISPHERE',' AUGER ELE
     &CTRON DIFFRACTION')
 529  FORMAT(///,14X,'TYPE OF CALCULATION : AZIMUTHAL LEED',' VARIATIONS
     &')
 530  FORMAT(///,11X,'TYPE OF CALCULATION : FULL HEMISPHERE',' LEED')
 531  FORMAT(///,18X,'TYPE OF CALCULATION : LEED ENERGY ','VARIATIONS')
 532  FORMAT(///,16X,'TYPE OF CALCULATION : POLAR LEED',' VARIATIONS')
 533  FORMAT(///,17X,'TYPE OF CALCULATION : EIGENVALUE',' ANALYSIS')
 534  FORMAT(///,22X,'THE AMPLITUDES WILL BE PRINTED SEPARATELY')
C
 701  FORMAT(6X,I1,1X,I3,2X,I4)
 702  FORMAT(6X,I1,1X,I3,3X,I3)
 703  FORMAT(15X,F8.3,3X,F8.3)
 713  FORMAT(6X,I1,1X,I3)
C
      RETURN
C
      END
