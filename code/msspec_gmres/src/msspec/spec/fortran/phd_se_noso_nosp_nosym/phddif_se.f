C
C=======================================================================
C
      SUBROUTINE PHDDIF_SE(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOK,NATC
     &LU,NFICHLEC,JFICH,NP)
C
C   This subroutine computes the PhD formula in the spin-independent case
C        from a non spin-orbit resolved initial core state LI.
C
C   Alternatively, it can compute the PhD amplitude for the APECS process.
C
C   The calculation is performed using a series expansion for the
C           expression of the scattering path operator
C
C                                        Last modified : 10 Jan 2016
C
      USE DIM_MOD
      USE PATH_SAMPLING_MOD
C
      USE ALGORITHM_MOD
      USE AMPLI_MOD
      USE APPROX_MOD
      USE COOR_MOD , NTCLU => NATCLU, NTP => NATYP
      USE C_RENORM_MOD
      USE DEBWAL_MOD
      USE DIRECT_MOD , RTHETA => RTHEXT
      USE EXTREM_MOD
      USE FIXSCAN_MOD
      USE INFILES_MOD
      USE INUNITS_MOD
      USE INIT_L_MOD
      USE INIT_J_MOD
      USE LIMAMA_MOD
      USE LINLBD_MOD
      USE MOYEN_MOD
      USE OUTFILES_MOD
      USE OUTUNITS_MOD
      USE PARCAL_MOD
      USE PATH_MOD
      USE PRINTP_MOD
      USE RENORM_MOD
      USE RESEAU_MOD
      USE SPIN_MOD
      USE TESTPA_MOD
      USE TESTPB_MOD
      USE TESTS_MOD
      USE TRANS_MOD
      USE TYPCAL_MOD
      USE TYPEM_MOD
      USE TYPEXP_MOD
      USE VALIN_MOD , PHLUM => PHILUM
      USE VALIN_AV_MOD
      USE VALFIN_MOD
C
      REAL NPATH1(0:NDIF_M),NOPA
      REAL LUM(3),AXE(3),EPS(3),DIRLUM(3),E_PH(NE_M)
C
      COMPLEX IC,ONEC,ZEROC,COEF,PW(0:NDIF_M),DELTA
      COMPLEX TLT(0:NT_M,4,NATM,NE_M),RHOMI
      COMPLEX TAU(LINMAX,LINFMAX,NATCLU_M)
      COMPLEX YLMR(0:NL_M,-NL_M:NL_M),MATRIX(3,2)
      COMPLEX YLME(0:NL_M,-NL_M:NL_M)
      COMPLEX R2,MLFLI(2,-LI_M:LI_M,3,2,3)
      COMPLEX SJDIR_1,SJDIR_2,SJDIF_1,SJDIF_2
      COMPLEX RHOK(NE_M,NATM,0:18,5,NSPIN2_M),RD
      COMPLEX SLJDIF,ATT_M,MLIL0(2,-LI_M:LI_M,6),SLF_1,SLF_2
      COMPLEX SL0DIF,SMJDIF
C
      DIMENSION VAL(NATCLU_M),NATYP(NATM),DIRPOL(3,2)
      DIMENSION EMET(3),R_L(9),COORD(3,NATCLU_M)
      DIMENSION R(NDIF_M),XR(NDIF_M),YR(NDIF_M),ZR(NDIF_M)
      DIMENSION JPOS(NDIF_M,3),JPA(NDIF_M)
C
C
C
      CHARACTER*7 STAT
      CHARACTER*13 OUTDATA1,OUTDATA2
C
C
      CHARACTER*24 OUTFILE
      CHARACTER*24 AMPFILE
C
      DATA PI,PIS180,CONV /3.141593,0.017453,0.512314/
      DATA FINSTRUC,CVECT,SMALL /0.007297,1.0,0.0001/
C
      CALL INIT_PATH_SAMPLING()
C
      ALGO1='SE'
      ALGO2='  '
      ALGO3='  '
      ALGO4='  '
C
      I_DIR=0
      NSET=1
      JEL=1
      OUTDATA1='CROSS-SECTION'
      IF(I_AMP.EQ.1) THEN
        I_MI=1
        OUTDATA2='MS AMPLITUDES'
      ELSE
        I_MI=0
      ENDIF
C
      IF(SPECTRO.EQ.'PHD') THEN
        IOUT=IUO2
        OUTFILE=OUTFILE2
        STAT='UNKNOWN'
        IF(I_MI.EQ.1) THEN
          IOUT2=IUSCR2+1
          N_DOT=1
          DO J_CHAR=1,24
            IF(OUTFILE(J_CHAR:J_CHAR).EQ.'.') GOTO 888
            N_DOT=N_DOT+1
          ENDDO
  888     CONTINUE
          AMPFILE=OUTFILE(1:N_DOT)//'amp'
          OPEN(UNIT=IOUT2, FILE=AMPFILE, STATUS=STAT)
        ENDIF
      ELSEIF(SPECTRO.EQ.'APC') THEN
        IOUT=IUSCR2+1
        OUTFILE='res/phot.amp'
        STAT='UNKNOWN'
      ENDIF
C
C  Position of the light when the analyzer is along the z axis :
C                  (X_LUM_Z,Y_LUM_Z,Z_LUM_Z)
C
      RTHLUM=THLUM*PIS180
      RPHLUM=PHLUM*PIS180
      X_LUM_Z=SIN(RTHLUM)*COS(RPHLUM)
      Y_LUM_Z=SIN(RTHLUM)*SIN(RPHLUM)
      Z_LUM_Z=COS(RTHLUM)
C
      IF(IMOD.EQ.0) THEN
C
C  The analyzer is rotated
C
        DIRLUM(1)=X_LUM_Z
        DIRLUM(2)=Y_LUM_Z
        DIRLUM(3)=Z_LUM_Z
      ELSE
C
C  The sample is rotated ---> light and analyzer rotated
C
        IF(I_EXT.EQ.0) THEN
          RTH0=THETA0*PIS180
          RPH0=PHI0*PIS180
          RTH=RTH0
          RPH=RPH0
C
C  R_L is the rotation matrix from 0z to (THETA0,PHI0) expressed as
C    a function of the Euler angles ALPHA=PHI0, BETA=THETA0, GAMMA=-PHI0
C    It is stored as (1 2 3)
C                    (4 5 6)
C                    (7 8 9)
C
          R_L(1)=COS(RTH0)*COS(RPH0)*COS(RPH0)+SIN(RPH0)*SIN(RPH0)
          R_L(2)=COS(RTH0)*SIN(RPH0)*COS(RPH0)-SIN(RPH0)*COS(RPH0)
          R_L(3)=SIN(RTH0)*COS(RPH0)
          R_L(4)=COS(RTH0)*SIN(RPH0)*COS(RPH0)-SIN(RPH0)*COS(RPH0)
          R_L(5)=COS(RTH0)*SIN(RPH0)*SIN(RPH0)+COS(RPH0)*COS(RPH0)
          R_L(6)=SIN(RTH0)*SIN(RPH0)
          R_L(7)=-SIN(RTH0)*COS(RPH0)
          R_L(8)=-SIN(RTH0)*SIN(RPH0)
          R_L(9)=COS(RTH0)
C
C  Position of the light when the analyzer is along (THETA0,PHI0) : LUM(3)
C
          LUM(1)=X_LUM_Z*R_L(1)+Y_LUM_Z*R_L(2)+Z_LUM_Z*R_L(3)
          LUM(2)=X_LUM_Z*R_L(4)+Y_LUM_Z*R_L(5)+Z_LUM_Z*R_L(6)
          LUM(3)=X_LUM_Z*R_L(7)+Y_LUM_Z*R_L(8)+Z_LUM_Z*R_L(9)
C
        ENDIF
      ENDIF
C
      IC=(0.,1.)
      ONEC=(1.,0.)
      ZEROC=(0.,0.)
      NSCAT=NATCLU-1
      ATTSE=1.
      ATTSJ=1.
      NPATH2(0)=1.
      NPATH(0)=1.
      NPMA(0)=1.
      NPMI(0)=1.
      ZSURF=VAL(1)
C
      IF((ISOM.EQ.0).OR.(JFICH.EQ.1)) THEN
        OPEN(UNIT=IOUT, FILE=OUTFILE, STATUS=STAT)
      ENDIF
C
C  Writing the headers in the output file
C
      CALL HEADERS(IOUT)
C
      IF((ISOM.EQ.0).OR.((ISOM.GT.0).AND.(JFICH.EQ.1))) THEN
        WRITE(IOUT,12) SPECTRO,OUTDATA1
        WRITE(IOUT,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE,IP
     &H_1,I_EXT
        IF(I_MI.EQ.1) THEN
          WRITE(IOUT2,12) SPECTRO,OUTDATA2
          WRITE(IOUT2,12) STEREO
          WRITE(IOUT2,19) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,I
     &E,IPH_1,I_EXT
          WRITE(IOUT2,20) PHI0,THETA0,PHI1,THETA1,NONVOL(1)
        ENDIF
      ENDIF
C
      IF(ISOM.EQ.0) THEN
        WRITE(IOUT,79) NPLAN,NEMET,NTHETA,NPHI,NE
        IF(I_MI.EQ.1) THEN
          WRITE(IOUT2,79) NPLAN,NEMET,NTHETA,NPHI,NE
        ENDIF
      ELSEIF((ISOM.NE.0).AND.(JFICH.EQ.1)) THEN
        WRITE(IOUT,11) NTHETA,NPHI,NE
        IF(I_MI.EQ.1) THEN
          WRITE(IOUT2,11) NTHETA,NPHI,NE
        ENDIF
      ENDIF
C
C   Construction of the linear index LAMBDA=(MU,NU)
C
      LAMBDA0=0
      DO N_O=0,NO
        NMX=N_O/2
        DO NU=0,NMX
          DO MU=-N_O,N_O
            NMU=2*NU+ABS(MU)
            IF(NMU.EQ.N_O) THEN
              LAMBDA0=LAMBDA0+1
              LBD(MU,NU)=LAMBDA0
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      LBDMAX=LAMBDA0
      IJK=0
C
C   Loop over the planes
C
      DO JPLAN=1,NPLAN
       Z=VAL(JPLAN)
       IF((IPHA.EQ.1).OR.(IPHA.EQ.2)) THEN
         DZZEM=ABS(Z-ZEM)
         IF(DZZEM.LT.SMALL) GOTO 10
         GOTO 1
       ENDIF
 10    CONTINUE
C
C   Loop over the different absorbers in a given plane
C
       DO JEMET=1,NEMET
        CALL EMETT(JEMET,IEMET,Z,SYM_AT,NATYP,EMET,NTYPEM,JNEM,*4)
        GO TO 2
   4    IF((ISORT1.EQ.0).AND.(IPRINT.GT.0)) THEN
          IF(I_TEST.NE.2) WRITE(IUO1,51) JPLAN,NTYPEM
        ENDIF
        GO TO 3
   2    IF((ABS(EMET(3)).GT.COUPUR).AND.(IBAS.EQ.1)) GOTO 5
        IF((ISORT1.EQ.0).AND.(IPRINT.GT.0)) THEN
          IF(I_TEST.NE.2) THEN
            WRITE(IUO1,52) JPLAN,EMET(1),EMET(2),EMET(3),NTYPEM
          ENDIF
        ENDIF
        IF(ISOM.EQ.1) NP=JPLAN
        ZSURFE=VAL(1)-EMET(3)
C
C   Loop over the energies
C
        DO JE=1,NE
         FMIN(0)=1.
         FMAX(0)=1.
         IF(NE.GT.1) THEN
           ECIN=E0+FLOAT(JE-1)*(EFIN-E0)/FLOAT(NE-1)
           E_PH(JE)=ELUM+FLOAT(JE-1)*(EFIN-E0)/FLOAT(NE-1)
         ELSEIF(NE.EQ.1) THEN
           ECIN=E0
           E_PH(JE)=ELUM
         ENDIF
         IF(I_TEST.NE.1) THEN
           CFM=8.*PI*E_PH(JE)*FINSTRUC
         ELSE
           CFM=1.
         ENDIF
         CALL LPM(ECIN,XLPM,*6)
         XLPM1=XLPM/A
         IF(IPRINT.GT.0) WRITE(IUO1,56) A,XLPM1
         IF((IPRINT.GT.0).AND.(IBAS.EQ.1)) THEN
           IF(I_TEST.NE.2) WRITE(IUO1,57) COUPUR
         ENDIF
         IF(ITL.EQ.0) THEN
           VK(JE)=SQRT(ECIN+VINT)*CONV*A*(1.,0.)
           VK2(JE)=CABS(VK(JE)*VK(JE))
         ENDIF
         GAMMA=1./(2.*XLPM1)
         IF(IPOTC.EQ.0) THEN
           VK(JE)=VK(JE)+IC*GAMMA
         ENDIF
         IF(I_TEST.NE.1) THEN
           VKR=REAL(VK(JE))
         ELSE
           VKR=1.
         ENDIF
         IF(I_MI.EQ.1) THEN
           WRITE(IOUT2,21) ECIN,VKR*CFM
         ENDIF
         IF((IDWSPH.EQ.1).AND.(ISPEED.EQ.1)) THEN
           IF(IDCM.GE.1) WRITE(IUO1,22)
           DO JAT=1,N_PROT
             IF(IDCM.EQ.0) THEN
               XK2UJ2=VK2(JE)*UJ2(JAT)
             ELSE
               XK2UJ2=VK2(JE)*UJ_SQ(JAT)
               WRITE(IUO1,23) JAT,UJ_SQ(JAT)*A*A
             ENDIF
             CALL DWSPH(JAT,JE,XK2UJ2,TLT,ISPEED)
             DO LAT=0,LMAX(JAT,JE)
               TL(LAT,1,JAT,JE)=TLT(LAT,1,JAT,JE)
             ENDDO
           ENDDO
         ENDIF
         IF(ABS(I_EXT).GE.1) THEN
           OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
           READ(IUI6,13) I_DIR,NSET,N_DUM1
           READ(IUI6,14) I_DUM1,N_DUM2,N_DUM3
         ENDIF
C
C   Initialization of TAU(INDJ,LINFMAX,JTYP)
C
         JATL=0
         DO JTYP=1,N_PROT
           NBTYP=NATYP(JTYP)
           LMJ=LMAX(JTYP,JE)
           DO JNUM=1,NBTYP
             JATL=JATL+1
             DO LF=LF1,LF2,ISTEP_LF
               ILF=LF*LF+LF+1
               DO MF=-LF,LF
                 INDF=ILF+MF
                 DO LJ=0,LMJ
                   ILJ=LJ*LJ+LJ+1
                   DO MJ=-LJ,LJ
                     INDJ=ILJ+MJ
                     TAU(INDJ,INDF,JATL)=ZEROC
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
         ENDDO
C
C   Storage of the coupling matrix elements MLFLI along the basis
C                       directions X,Y ET Z
C
C   These basis directions refer to the polarization if IDICHR = 0
C                 but to the light when IDICHR = 1
C
C                      JBASE = 1 : X
C                      JBASE = 2 : Y
C                      JBASE = 3 : Z
C
         DO MI=-LI,LI
           DO LF=LF1,LF2,ISTEP_LF
             LR=1+(1+LF-LI)/2
             DELTA=DLT(JE,NTYPEM,NNL,LR)
             RD=RHOK(JE,NTYPEM,NNL,LR,1)
             DO MF=-LF,LF
               IF((MF.LT.MI-1).OR.(MF.GT.MI+1)) GOTO 333
               IF((INITL.EQ.0).AND.(MF.NE.MI)) GOTO 333
               MR=2+MF-MI
               CALL COUMAT(ITL,MI,LF,MF,DELTA,RD,MATRIX)
               DO JBASE=1,3
                 MLFLI(1,MI,MR,LR,JBASE)=MATRIX(JBASE,1)
                 IF(IDICHR.GE.1) THEN
                   MLFLI(2,MI,MR,LR,JBASE)=MATRIX(JBASE,2)
                 ENDIF
               ENDDO
 333           CONTINUE
             ENDDO
           ENDDO
         ENDDO
C
C    Calculation of the scattering path operator TAU
C
         IF(I_TEST.EQ.2) GOTO 666
         PW(0)=ONEC
         PW(1)=ONEC
         ND=0
         TH01=0.
         PHI01=0.
         RHO01=ZEROC
         THMI=0.
         PHMI=0.
         RHOMI=ZEROC
         JATLEM=JNEM
         IF(NTYPEM.GT.1) THEN
           DO JAEM=NTYPEM-1,1,-1
             JATLEM=JATLEM+NATYP(JAEM)
           ENDDO
         ENDIF
         DO JD=1,NDIF
           NPATH2(JD)=0.
           NPATH(JD)=0.
           IT(JD)=0
           IN(JD)=0
           FMIN(JD)=1.E+20
           FMAX(JD)=0.
         ENDDO
         NTHOF=0
C
C     Calculation of the maximal intensity for the paths of order NCUT
C      (plane waves). This will be taken as a reference for the IPW filter.
C
         IF(IPW.EQ.1) THEN
           NDIFOLD=NDIF
           NOOLD=NO
           ISPHEROLD=ISPHER
           NDIF=NCUT
           NO=0
           ISPHER=0
           IREF=1
           IPW=0
           IJ=0
           DIJ=0.
           FREF=0.
           CALL FINDPATHS(ND,NTYPEM,JATLEM,I_CP,R,XR,YR,ZR,RHOMI,THMI,PH
     &MI,ZSURF,JPOS,PW,JE,FREF,DIJ,TAU)
           NDIF=NDIFOLD
           NO=NOOLD
           ISPHER=ISPHEROLD
           PW(0)=ONEC
           PW(1)=ONEC
           IPW=1
           ND=0
           TH01=0.
           PHI01=0.
           RHO01=ZEROC
           THMI=0.
           PHMI=0.
           RHOMI=ZEROC
           JATLEM=JNEM
           IF(NTYPEM.GT.1) THEN
             DO JAEM=NTYPEM-1,1,-1
               JATLEM=JATLEM+NATYP(JAEM)
             ENDDO
           ENDIF
           DO JD=1,NDIF
             NPATH2(JD)=0.
             NPATH(JD)=0.
             IT(JD)=0
             IN(JD)=0
             FMIN(JD)=1.E+20
             FMAX(JD)=0.
           ENDDO
           NTHOF=0
C
C   New initialization of TAU(INDJ,INDF,JATL) after the PW calculation
C
           JATL=0
           DO JTYP=1,N_PROT
             NBTYP=NATYP(JTYP)
             LMJ=LMAX(JTYP,JE)
             DO JNUM=1,NBTYP
               JATL=JATL+1
               DO LF=LF1,LF2,ISTEP_LF
                 ILF=LF*LF+LF+1
                 DO MF=-LF,LF
                   INDF=ILF+MF
                   DO LJ=0,LMJ
                     ILJ=LJ*LJ+LJ+1
                     DO MJ=-LJ,LJ
                       INDJ=ILJ+MJ
                       TAU(INDJ,INDF,JATL)=ZEROC
                     ENDDO
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
         ENDIF
C
C   Generation and print-out of the paths
C
         IF (NPATHP.GT.0) THEN
           DO JP=1,NPATHP-1
             FMN(JP)=0.
             PATH(JP)=0.
             JON(JP)=0
           ENDDO
           FMN(NPATHP)=-1.
           PATH(NPATHP)=0.
           JON(NPATHP)=0
         ENDIF
         IREF=0
         IJ=1
         IF(IPRINT.EQ.3) THEN
           OPEN(UNIT=IUSCR, STATUS='SCRATCH')
         ENDIF
         DIJ=0.
         CALL FINDPATHS(ND,NTYPEM,JATLEM,I_CP,R,XR,YR,ZR,RHOMI,THMI,PHMI
     &,ZSURF,JPOS,PW,JE,FREF,DIJ,TAU)
         IF(NPATHP.EQ.0) GOTO 15
         IF(NSCAT.GT.1) THEN
           XPATOT=REAL((DFLOAT(NSCAT)**DFLOAT(NDIF+1) -1.D0)/DFLOAT(NSCA
     &T-1))
         ELSE
           XPATOT=FLOAT(NDIF+1)
         ENDIF
         IF(XPATOT.LT.2.14748E+09) THEN
           NPATOT=INT(XPATOT)
           IF(NPATOT.LT.NPATHP) NPATHP=NPATOT-1
         ENDIF
         WRITE(IUO1,84) NPATHP
         WRITE(IUO1,81)
         DO JPT=1,NPATHP
           IF(PATH(NPATHP).GT.2.14E+09) THEN
CST           WRITE(IUO1,82) JPT,JON(JPT),PATH(JPT),FMN(JPT),DMN(JPT),JNE
CST  &M,(JPON(JPT,KD),KD=1,JON(JPT))
             PRINT *, JPT,JON(JPT),PATH(JPT),FMN(JPT),DMN(JPT),JNE
     &M,(JPON(JPT,KD),KD=1,JON(JPT))
           ELSE
CST         WRITE(IUO1,83) JPT,JON(JPT),INT(PATH(JPT)),FMN(JPT),DMN(JPT
CST &),JNEM,(JPON(JPT,KD),KD=1,JON(JPT))
             PRINT *, JPT,JON(JPT),INT(PATH(JPT)),FMN(JPT),DMN(JPT
     &),JNEM,(JPON(JPT,KD),KD=1,JON(JPT))
           ENDIF
         ENDDO
         IF(IPRINT.EQ.3) THEN
           IF(XPATOT.GT.2.14748E+09) GOTO 172
           WRITE(IUO1,85)
           WRITE(IUO1,71)
           NPATOT=INT(XPATOT)
           DO JOP=0,NDIF
             IF(JOP.EQ.0) THEN
               XINT0=FMAX(0)
               DIST0=0.
               WRITE(IUO1,70) JOP,JOP+1,XINT0,DIST0,JNEM
               GOTO 75
             ENDIF
             WRITE(IUO1,77)
             DO JLINE=1,NPATOT-1
               READ(IUSCR,69,ERR=75,END=75) JOPA,NOPA,XMAX,DIST0,(JPA(KD
     &),KD=1,JOPA)
               IF(JOPA.EQ.JOP) THEN
                 IF(NOPA.GT.2.14E+09) THEN
                  WRITE(IUO1,76) JOPA,NOPA,XMAX,DIST0,JNEM,(JPA(KD),KD=1
     &,JOPA)
                 ELSE
                  WRITE(IUO1,70) JOPA,INT(NOPA),XMAX,DIST0,JNEM,(JPA(KD)
     &,KD=1,JOPA)
                 ENDIF
               ENDIF
             ENDDO
             IF(JOP.EQ.NDIF) WRITE(IUO1,80)
  75         REWIND IUSCR
           ENDDO
           GOTO 73
  172      WRITE(IUO1,74)
           CLOSE(IUSCR,STATUS='DELETE')
  73     ENDIF
         DO JD=0,NDIF
           NPATH1(JD)=REAL(DFLOAT(NSCAT)**DFLOAT(JD))
           IF(NPATH1(JD).GT.2.14E+09) THEN
             IF(FMIN(JD).EQ.0.1E+21) FMIN(JD)=0.
             WRITE(IUO1,53) JD,NPATH1(JD),NPATH2(JD),FMIN(JD),NPMI(JD),F
     &MAX(JD),NPMA(JD)
             IF((IPW.EQ.1).AND.(JD.GT.NCUT)) WRITE(IUO1,68) FREF*PCTINT
           ELSE
             IF(FMIN(JD).EQ.0.1E+21) FMIN(JD)=0.
             WRITE(IUO1,58) JD,INT(NPATH1(JD)+0.1),INT(NPATH2(JD)+0.1),F
     &MIN(JD),INT(NPMI(JD)+0.1),FMAX(JD),INT(NPMA(JD)+0.1)
             IF((IPW.EQ.1).AND.(JD.GT.NCUT)) WRITE(IUO1,68) FREF*PCTINT
           ENDIF
         ENDDO
 666     CONTINUE
C
C    Calculation of the Photoelectron Diffraction formula
C
C
C    Loop over the 'fixed' angle
C
  15     DO J_FIXED=1,N_FIXED
          IF(N_FIXED.GT.1) THEN
            IF(I_EXT.EQ.0) THEN
              FIX_STEP=(FIX1-FIX0)/FLOAT(N_FIXED-1)
              XINCRF=FLOAT(J_FIXED-1)*FIX_STEP
            ELSE
              XINCRF=0.
            ENDIF
          ELSEIF(N_FIXED.EQ.1) THEN
            XINCRF=0.
          ENDIF
          IF(ABS(I_EXT).GE.1) THEN
            READ(IUI6,86) JSET,JLINE,THD,PHD
            IF(I_EXT.EQ.-1) BACKSPACE IUI6
            THETA0=THD
            PHI0=PHD
          ENDIF
          IF(IPH_1.EQ.1) THEN
             IF(I_EXT.EQ.0) THEN
               DPHI=PHI0+XINCRF
             ELSE
               DPHI=PHD
             ENDIF
             RPHI=DPHI*PIS180
             IF(IPRINT.GT.0) WRITE(IUO1,66) DPHI
          ELSE
             ISAUT=0
             IF(I_EXT.EQ.0) THEN
               DTHETA=THETA0+XINCRF
             ELSE
               DTHETA=THD
             ENDIF
             RTHETA=DTHETA*PIS180
             IF(ABS(DTHETA).GT.90.) ISAUT=ISAUT+1
             IF(I_EXT.GE.1) ISAUT=0
             IF(I_TEST.EQ.2)  ISAUT=0
             IF(ISAUT.GT.0) GOTO 8
             IF(IPRINT.GT.0) WRITE(IUO1,65) DTHETA
             IF((IPRINT.GT.0).AND.(I_TEST.NE.2)) WRITE(IUO1,59)
             IF((IPRINT.EQ.1).AND.(I_TEST.NE.2)) WRITE(IUO1,60)
C
C  THETA-dependent number of PHI points for stereographic
C     representation (to obtain a uniform sampling density).
C     (Courtesy of J. Osterwalder - University of Zurich)
C
             IF(STEREO.EQ.'YES') THEN
               N_SCAN=INT((SCAN1-SCAN0)*SIN(RTHETA)/FIX_STEP+SMALL)+1
             ENDIF
C
          ENDIF
          IF((N_FIXED.GT.1).AND.(IMOD.EQ.1)) THEN
C
C  When there are several sets of scans (N_FIXED > 1),
C   the initial position LUM of the light is recalculated
C      for each initial position (RTH,RPH) of the analyzer
C
            IF(IPH_1.EQ.1) THEN
                RTH=THETA0*PIS180
                RPH=RPHI
            ELSE
              RTH=RTHETA
              RPH=PHI0*PIS180
            ENDIF
C
            R_L(1)=COS(RTH)*COS(RPH)
            R_L(2)=-SIN(RPH)
            R_L(3)=SIN(RTH)*COS(RPH)
            R_L(4)=COS(RTH)*SIN(RPH)
            R_L(5)=COS(RPH)
            R_L(6)=SIN(RTH)*SIN(RPH)
            R_L(7)=-SIN(RTH)
            R_L(8)=0.
            R_L(9)=COS(RTH)
C
            LUM(1)=X_LUM_Z*R_L(1)+Y_LUM_Z*R_L(2)+Z_LUM_Z*R_L(3)
            LUM(2)=X_LUM_Z*R_L(4)+Y_LUM_Z*R_L(5)+Z_LUM_Z*R_L(6)
            LUM(3)=X_LUM_Z*R_L(7)+Y_LUM_Z*R_L(8)+Z_LUM_Z*R_L(9)
          ENDIF
C
C    Loop over the scanned angle
C
          DO J_SCAN=1,N_SCAN
           IF(N_SCAN.GT.1) THEN
             XINCRS=FLOAT(J_SCAN-1)*(SCAN1-SCAN0)/FLOAT(N_SCAN-1)
           ELSEIF(N_SCAN.EQ.1) THEN
             XINCRS=0.
           ENDIF
           IF(I_EXT.EQ.-1) THEN
             READ(IUI6,86) JSET,JLINE,THD,PHD
             BACKSPACE IUI6
           ENDIF
           IF(IPH_1.EQ.1) THEN
             ISAUT=0
             IF(I_EXT.EQ.0) THEN
               DTHETA=THETA0+XINCRS
             ELSE
               DTHETA=THD
             ENDIF
             RTHETA=DTHETA*PIS180
             IF(ABS(DTHETA).GT.90.) ISAUT=ISAUT+1
             IF(I_EXT.GE.1) ISAUT=0
             IF(I_TEST.EQ.2)  ISAUT=0
             IF(ISAUT.GT.0) GOTO 8
             IF(IPRINT.GT.0) WRITE(IUO1,65) DTHETA
             IF((IPRINT.GT.0).AND.(I_TEST.NE.2)) WRITE(IUO1,59)
             IF((IPRINT.EQ.1).AND.(I_TEST.NE.2)) WRITE(IUO1,60)
           ELSE
             IF(I_EXT.EQ.0) THEN
               DPHI=PHI0+XINCRS
             ELSE
               DPHI=PHD
             ENDIF
             RPHI=DPHI*PIS180
             IF(IPRINT.GT.0) WRITE(IUO1,66) DPHI
           ENDIF
C
C  Loop over the sets of directions to average over (for gaussian average)
C
C
           SSETDIR_1=0.
           SSETDIF_1=0.
           SSETDIR_2=0.
           SSETDIF_2=0.
C
           SSET2DIR_1=0.
           SSET2DIF_1=0.
           SSET2DIR_2=0.
           SSET2DIF_2=0.
C
           IF(I_EXT.EQ.-1) THEN
             JREF=INT(NSET)/2+1
           ELSE
             JREF=1
           ENDIF
C
           DO J_SET=1,NSET
           IF(I_EXT.EQ.-1) THEN
             READ(IUI6,86) JSET,JLINE,THD,PHD,W
             DTHETA=THD
             DPHI=PHD
             RTHETA=DTHETA*PIS180
             RPHI=DPHI*PIS180
C
C  Here, there are several sets of scans (NSET > 1), so
C   the initial position LUM of the light must be
C   recalculated for each initial position of the analyzer
C
             RTH=TH_0(J_SET)*PIS180
             RPH=PH_0(J_SET)*PIS180
C
             IF(IMOD.EQ.1) THEN
               R_L(1)=COS(RTH)*COS(RPH)
               R_L(2)=-SIN(RPH)
               R_L(3)=SIN(RTH)*COS(RPH)
               R_L(4)=COS(RTH)*SIN(RPH)
               R_L(5)=COS(RPH)
               R_L(6)=SIN(RTH)*SIN(RPH)
               R_L(7)=-SIN(RTH)
               R_L(8)=0.
               R_L(9)=COS(RTH)
C
               LUM(1)=X_LUM_Z*R_L(1)+Y_LUM_Z*R_L(2)+Z_LUM_Z*R_L(3)
               LUM(2)=X_LUM_Z*R_L(4)+Y_LUM_Z*R_L(5)+Z_LUM_Z*R_L(6)
               LUM(3)=X_LUM_Z*R_L(7)+Y_LUM_Z*R_L(8)+Z_LUM_Z*R_L(9)
C
             ENDIF
           ELSE
             W=1.
           ENDIF
C
           IF(I_EXT.EQ.-1) PRINT 89
C
           CALL DIRAN(VINT,ECIN,JEL)
C
           IF(J_SET.EQ.JREF) THEN
             DTHETAP=DTHETA
             DPHIP=DPHI
           ENDIF
C
           IF(I_EXT.EQ.-1) THEN
             WRITE(IUO1,88) DTHETA,DPHI
           ENDIF
C
C     ..........          Case IMOD=1 only          ..........
C
C  Calculation of the position of the light when the analyzer is at
C   (THETA,PHI). DIRLUM is the direction of the light and its initial
C   value (at (THETA0,PHI0)) is LUM. AXE is the direction of the theta
C   rotation axis and EPS is defined so that (AXE,DIRLUM,EPS) is a
C   direct orthonormal basis. The transform of a vector R by a rotation
C   of OMEGA about AXE is then given by
C
C     R' = R COS(OMEGA) + (AXE.R)(1-COS(OMEGA)) AXE + (AXE^R) SIN(OMEGA)
C
C   Here, DIRANA is the internal direction of the analyzer and ANADIR
C                      its external position
C
C   Note that when the initial position of the analyzer is (RTH,RPH)
C    which coincides with (RTH0,RPH0) only for the first fixed angle
C
           IF(IMOD.EQ.1) THEN
             IF(ITHETA.EQ.1) THEN
               AXE(1)=-SIN(RPH)
               AXE(2)=COS(RPH)
               AXE(3)=0.
               RANGLE=RTHETA-RTH
             ELSEIF(IPHI.EQ.1) THEN
               AXE(1)=0.
               AXE(2)=0.
               AXE(3)=1.
               RANGLE=RPHI-RPH
             ENDIF
             CALL PRVECT(AXE,LUM,EPS,CVECT)
             PRS=PRSCAL(AXE,LUM)
             IF(J_SCAN.EQ.1) THEN
               DIRLUM(1)=LUM(1)
               DIRLUM(2)=LUM(2)
               DIRLUM(3)=LUM(3)
             ELSE
               DIRLUM(1)=LUM(1)*COS(RANGLE)+PRS*(1.-COS(RANGLE))*AXE(1)+
     &SIN(RANGLE)*EPS(1)
               DIRLUM(2)=LUM(2)*COS(RANGLE)+PRS*(1.-COS(RANGLE))*AXE(2)+
     &SIN(RANGLE)*EPS(2)
               DIRLUM(3)=LUM(3)*COS(RANGLE)+PRS*(1.-COS(RANGLE))*AXE(3)+
     &SIN(RANGLE)*EPS(3)
             ENDIF
           ENDIF
           IF(DIRLUM(3).GT.1.) DIRLUM(3)=1.
           IF(DIRLUM(3).LT.-1.) DIRLUM(3)=-1.
           THETALUM=ACOS(DIRLUM(3))
           IF(I_TEST.EQ.2) THETALUM=-THETALUM
           COEF=DIRLUM(1)+IC*DIRLUM(2)
           CALL ARCSIN(COEF,DIRLUM(3),PHILUM)
           ANALUM=ANADIR(1,1)*DIRLUM(1) + ANADIR(2,1)*DIRLUM(2) +ANADIR(
     &3,1)*DIRLUM(3)
C
           SEPSDIR_1=0.
           SEPSDIF_1=0.
           SEPSDIR_2=0.
           SEPSDIF_2=0.
C
C    Loop over the directions of polarization
C
           DO JEPS=1,NEPS
            IF((JEPS.EQ.1).AND.(IPOL.GE.0)) THEN
              DIRPOL(1,JEPS)=COS(THETALUM)*COS(PHILUM)
              DIRPOL(2,JEPS)=COS(THETALUM)*SIN(PHILUM)
              DIRPOL(3,JEPS)=-SIN(THETALUM)
            ELSE
              DIRPOL(1,JEPS)=-SIN(PHILUM)
              DIRPOL(2,JEPS)=COS(PHILUM)
              DIRPOL(3,JEPS)=0.
            ENDIF
            IF(ABS(IPOL).EQ.1) THEN
              IF(IPRINT.GT.0) THEN
                WRITE(IUO1,61) (DIRANA(J,1),J=1,3),(DIRLUM(K),K=1,3),   
     &    (DIRPOL(K,1),K=1,3),ANALUM
              ENDIF
            ELSE
              IF((JEPS.EQ.1).AND.(IPRINT.GT.0)) THEN
                WRITE(IUO1,63) (DIRANA(J,1),J=1,3),(DIRLUM(K),K=1,3),ANA
     &LUM
              ENDIF
            ENDIF
            IF((JEPS.EQ.1).AND.(I_EXT.EQ.-1)) PRINT 89
C
C   Calculation of the coupling matrix MLIL0
C
            DO MI=-LI,LI
              DO LF=LF1,LF2,ISTEP_LF
                LR=1+(1+LF-LI)/2
                LRR=3*(LR-1)
                DO MF=-LF,LF
                  MR=2+MF-MI
                  IF((MF.LT.MI-1).OR.(MF.GT.MI+1)) GOTO 777
                  IF((INITL.EQ.0).AND.(MF.NE.MI)) GOTO 777
                  LMR=LRR+MR
                  IF(IDICHR.EQ.0) THEN
                    IF(I_TEST.NE.1) THEN
                      MLIL0(1,MI,LMR)=MLFLI(1,MI,MR,LR,1)*DIRPOL(1,JEPS)
     & +MLFLI(1,MI,MR,LR,2)*DIRPOL(2,JEPS) +MLFLI(1,MI,MR,LR,3)*DIRPOL(3
     &,JEPS)
                    ELSE
                      MLIL0(1,MI,LMR)=ONEC
                    ENDIF
                  ELSEIF(IDICHR.GE.1) THEN
                    IF(I_TEST.NE.1) THEN
                      MLIL0(1,MI,LMR)=MLFLI(1,MI,MR,LR,1)*DIRLUM(1) +MLF
     &LI(1,MI,MR,LR,2)*DIRLUM(2) +MLFLI(1,MI,MR,LR,3)*DIRLUM(3)
                      MLIL0(2,MI,LMR)=MLFLI(2,MI,MR,LR,1)*DIRLUM(1) +MLF
     &LI(2,MI,MR,LR,2)*DIRLUM(2) +MLFLI(2,MI,MR,LR,3)*DIRLUM(3)
                    ELSE
                      MLIL0(1,MI,LMR)=ONEC
                    ENDIF
                  ENDIF
  777             CONTINUE
                ENDDO
              ENDDO
            ENDDO
C
            SRDIF_1=0.
            SRDIR_1=0.
            SRDIF_2=0.
            SRDIR_2=0.
C
C    Loop over the different directions of the analyzer contained in a cone
C
            DO JDIR=1,NDIR
             IF(IATTS.EQ.1) THEN
               ATTSE=EXP(-ZSURFE*GAMMA/DIRANA(3,JDIR))
             ENDIF
C
             SMIDIR_1=0.
             SMIDIF_1=0.
             SMIDIR_2=0.
             SMIDIF_2=0.
C
C    Loop over the equiprobable azimuthal quantum numbers MI corresponding
C                          to the initial state LI
C
             LME=LMAX(1,JE)
             CALL HARSPH(NL_M,THETAR(JDIR),PHIR(JDIR),YLME,LME)
             DO MI=-LI,LI
              SJDIR_1=ZEROC
              SJDIF_1=ZEROC
              SJDIR_2=ZEROC
              SJDIF_2=ZEROC
C
C    Calculation of the direct emission (used a a reference for the
C      output), which is not contained in the calculation of TAU
C
              DO LF=LF1,LF2,ISTEP_LF
                LR=1+(1+LF-LI)/2
                LRR=3*(LR-1)
                ILF=LF*LF+LF+1
                IF(ISPEED.EQ.1) THEN
                  R2=TL(LF,1,1,JE)
                ELSE
                  R2=TLT(LF,1,1,JE)
                ENDIF
                IF(I_REN.GE.1) R2=R2*C_REN(0)
                DO MF=-LF,LF
                  MR=2+MF-MI
                  LMR=LRR+MR
                  INDF=ILF+MF
                  IF((MF.LT.MI-1).OR.(MF.GT.MI+1)) GOTO 444
                  IF((INITL.EQ.0).AND.(MF.NE.MI)) GOTO 444
                  SJDIR_1=SJDIR_1+YLME(LF,MF)*ATTSE*MLIL0(1,MI,LMR)*R2
                  IF(IDICHR.GE.1) THEN
                    SJDIR_2=SJDIR_2+YLME(LF,MF)*ATTSE*MLIL0(2,MI,LMR)*R2
     &
                  ENDIF
C
C    Contribution of the absorber to TAU (initialization of SJDIF)
C
                  IF(I_TEST.EQ.2) GOTO 444
                  SL0DIF=ZEROC
                  DO L0=0,LME
                    IL0=L0*L0+L0+1
                    SL0DIF=SL0DIF+YLME(L0,0)*TAU(IL0,INDF,1)
                    DO M0=1,L0
                     IND01=IL0+M0
                     IND02=IL0-M0
                     SL0DIF=SL0DIF+(YLME(L0,M0)*TAU(IND01,INDF,1)+YLME(L
     &0,-M0)*TAU(IND02,INDF,1))
                    ENDDO
                  ENDDO
                  SJDIF_1=SJDIF_1+SL0DIF*MLIL0(1,MI,LMR)
                  IF(IDICHR.GE.1) THEN
                    SJDIF_2=SJDIF_2+SL0DIF*MLIL0(2,MI,LMR)
                  ENDIF
 444              CONTINUE
                ENDDO
              ENDDO
              SJDIF_1=SJDIF_1*ATTSE
              IF(IDICHR.GE.1) THEN
                SJDIF_2=SJDIF_2*ATTSE
              ENDIF
C
C    Loop over the last atom J encountered by the photoelectron
C                  before escaping the solid
C
              IF(I_TEST.EQ.2) GOTO 111
              DO JTYP=2,N_PROT
               NBTYP=NATYP(JTYP)
               LMJ=LMAX(JTYP,JE)
               DO JNUM=1,NBTYP
                JATL=NCORR(JNUM,JTYP)
                XOJ=SYM_AT(1,JATL)-EMET(1)
                YOJ=SYM_AT(2,JATL)-EMET(2)
                ZOJ=SYM_AT(3,JATL)-EMET(3)
                ROJ=SQRT(XOJ*XOJ+YOJ*YOJ+ZOJ*ZOJ)
                ZSURFJ=VAL(1)-SYM_AT(3,JATL)
                CALL HARSPH(NL_M,THETAR(JDIR),PHIR(JDIR),YLMR,LMJ)
                IF(IATTS.EQ.1) THEN
                  ATTSJ=EXP(-ZSURFJ*GAMMA/DIRANA(3,JDIR))
                ENDIF
                CSTHJR=(XOJ*DIRANA(1,JDIR)+YOJ*DIRANA(2,JDIR)+ZOJ*DIRANA
     &(3,JDIR))/ROJ
                IF((IDWSPH.EQ.1).AND.(ISPEED.EQ.1)) GOTO 78
                CTROIS1=ZOJ/ROJ
                IF(CTROIS1.GT.1.)  THEN
                  CTROIS1=1.
                ELSEIF(CTROIS1.LT.-1.) THEN
                  CTROIS1=-1.
                ENDIF
                IF(IDCM.GE.1) THEN
                  UJ2(JTYP)=UJ_SQ(JTYP)
                ENDIF
                IF(ABS(ZSURFJ).LE.SMALL) THEN
                  IF(ABS(CSTHJR-1.).GT.SMALL) THEN
                    CSKZ2J=(DIRANA(3,JDIR)-CTROIS1)*(DIRANA(3,JDIR)-CTRO
     &IS1)/(2.-2.*CSTHJR)
                  ELSE
                    CSKZ2J=1.
                  ENDIF
                  UJJ=UJ2(JTYP)*(1.+CSKZ2J*(RSJ-1.))
                ELSE
                  UJJ=UJ2(JTYP)
                ENDIF
                IF((ISPEED.EQ.0).AND.(IDWSPH.EQ.1)) THEN
                  XK2UJ2=VK2(JE)*UJJ
                  CALL DWSPH(JTYP,JE,XK2UJ2,TLT,ISPEED)
                ENDIF
 78             IF(IDWSPH.EQ.1) THEN
                  DWTER=1.
                ELSE
                  DWTER=EXP(-VK2(JE)*UJJ*(1.-CSTHJR))
                ENDIF
                IF(JATL.EQ.JATLEM) THEN
                  ATT_M=ATTSE*DWTER
                ELSE
                  ATT_M=ATTSJ*DWTER*CEXP(-IC*VK(JE)*ROJ*CSTHJR)
                ENDIF
C
                SLF_1=ZEROC
                SLF_2=ZEROC
                DO LF=LF1,LF2,ISTEP_LF
                 LR=1+(1+LF-LI)/2
                 LRR=3*(LR-1)
                 ILF=LF*LF+LF+1
                 DO MF=-LF,LF
                  MR=2+MF-MI
                  INDF=ILF+MF
                  IF((MF.LT.MI-1).OR.(MF.GT.MI+1)) GOTO 555
                  IF((INITL.EQ.0).AND.(MF.NE.MI)) GOTO 555
                  LMR=LRR+MR
                  SLJDIF=ZEROC
                  DO LJ=0,LMJ
                   ILJ=LJ*LJ+LJ+1
                   SMJDIF=YLMR(LJ,0)*TAU(ILJ,INDF,JATL)
                   IF(LJ.GT.0) THEN
                     DO MJ=1,LJ
                      INDJ1=ILJ+MJ
                      INDJ2=ILJ-MJ
                      SMJDIF=SMJDIF+(YLMR(LJ,MJ)*TAU(INDJ1,INDF,JATL)+YL
     &MR(LJ,-MJ)*TAU(INDJ2,INDF,JATL))
                     ENDDO
                   ENDIF
                   SLJDIF=SLJDIF+SMJDIF
                  ENDDO
                  SLF_1=SLF_1+SLJDIF*MLIL0(1,MI,LMR)
                  IF(IDICHR.GE.1) THEN
                    SLF_2=SLF_2+SLJDIF*MLIL0(2,MI,LMR)
                  ENDIF
 555              CONTINUE
                 ENDDO
                ENDDO
                SJDIF_1=SJDIF_1+SLF_1*ATT_M
                IF(IDICHR.GE.1) THEN
                  SJDIF_2=SJDIF_2+SLF_2*ATT_M
                ENDIF
C
C     End of the loops over the last atom J
C
               ENDDO
              ENDDO
C
C     Writing the amplitudes in file IOUT for APECS, or
C       in file IOUT2 for PhD (orientated orbitals' case)
C
 111          IF(SPECTRO.EQ.'APC') THEN
                WRITE(IOUT,87) JFICH,JPLAN,JEMET,JE,J_FIXED,J_SCAN,JEPS,
     &JDIR,MI,SJDIR_1,SJDIR_1+SJDIF_1
                IF(IDICHR.GE.1) THEN
                  WRITE(IOUT,87) JFICH,JPLAN,JEMET,JE,J_FIXED,J_SCAN,JEP
     &S,JDIR,MI,SJDIR_2,SJDIR_2+SJDIF_2
                ENDIF
              ELSE
                IF(I_MI.EQ.1) THEN
                  WRITE(IOUT2,87) JFICH,JPLAN,JEMET,JE,J_FIXED,J_SCAN,JE
     &PS,JDIR,MI,SJDIR_1,SJDIR_1+SJDIF_1
                  IF(IDICHR.GE.1) THEN
                    WRITE(IOUT2,87) JFICH,JPLAN,JEMET,JE,J_FIXED,J_SCAN,
     &JEPS,JDIR,MI,SJDIR_2,SJDIR_2+SJDIF_2
                  ENDIF
                ENDIF
C
C     Computing the square modulus
C
                SMIDIF_1=SMIDIF_1+CABS(SJDIR_1+SJDIF_1)*CABS(SJDIR_1+SJD
     &IF_1)
                SMIDIR_1=SMIDIR_1+CABS(SJDIR_1)*CABS(SJDIR_1)
                IF(IDICHR.GE.1) THEN
                  SMIDIF_2=SMIDIF_2+CABS(SJDIR_2+SJDIF_2)*CABS(SJDIR_2+S
     &JDIF_2)
                  SMIDIR_2=SMIDIR_2+CABS(SJDIR_2)*CABS(SJDIR_2)
                ENDIF
              ENDIF
C
C     End of the loop over MI
C
             ENDDO
C
             IF(SPECTRO.EQ.'APC') GOTO 220
             SRDIR_1=SRDIR_1+SMIDIR_1
             SRDIF_1=SRDIF_1+SMIDIF_1
             IF(IDICHR.GE.1) THEN
               SRDIR_2=SRDIR_2+SMIDIR_2
               SRDIF_2=SRDIF_2+SMIDIF_2
             ENDIF
 220         CONTINUE
C
C     End of the loop on the directions of the analyzer
C
            ENDDO
C
            IF(SPECTRO.EQ.'APC') GOTO 221
            SEPSDIF_1=SEPSDIF_1+SRDIF_1*VKR*CFM/NDIR
            SEPSDIR_1=SEPSDIR_1+SRDIR_1*VKR*CFM/NDIR
            IF(IDICHR.GE.1) THEN
              SEPSDIF_2=SEPSDIF_2+SRDIF_2*VKR*CFM/NDIR
              SEPSDIR_2=SEPSDIR_2+SRDIR_2*VKR*CFM/NDIR
            ENDIF
 221        CONTINUE
C
C     End of the loop on the polarization
C
           ENDDO
C
           SSETDIR_1=SSETDIR_1+SEPSDIR_1*W
           SSETDIF_1=SSETDIF_1+SEPSDIF_1*W
           IF(ICHKDIR.EQ.2) THEN
             IF(JSET.EQ.JREF) THEN
               SSET2DIR_1=SEPSDIR_1
               SSET2DIF_1=SEPSDIF_1
             ENDIF
           ENDIF
           IF(IDICHR.GE.1) THEN
             SSETDIR_2=SSETDIR_2+SEPSDIR_2*W
             SSETDIF_2=SSETDIF_2+SEPSDIF_2*W
             IF(ICHKDIR.EQ.2) THEN
               IF(JSET.EQ.JREF) THEN
                 SSET2DIR_2=SEPSDIR_2
                 SSET2DIF_2=SEPSDIF_2
               ENDIF
             ENDIF
           ENDIF
C
C     End of the loop on the set averaging
C
           ENDDO
C
           IF(SPECTRO.EQ.'APC') GOTO 222
           IF(IDICHR.EQ.0) THEN
             IF(ISOM.EQ.2) THEN
              WRITE(IOUT,67) JPLAN,JFICH,DTHETAP,DPHIP,ECIN,SSETDIR_1,SS
     &ETDIF_1
              IF(ICHKDIR.EQ.2) THEN
               WRITE(IOUT,67) JPLAN,JFICH,DTHETAP,DPHIP,ECIN,SSET2DIR_1,
     &SSET2DIF_1
              ENDIF
             ELSE
              WRITE(IOUT,67) JPLAN,JEMET,DTHETAP,DPHIP,ECIN,SSETDIR_1,SS
     &ETDIF_1
              IF(ICHKDIR.EQ.2) THEN
               WRITE(IOUT,67) JPLAN,JEMET,DTHETAP,DPHIP,ECIN,SSET2DIR_1,
     &SSET2DIF_1
              ENDIF
             ENDIF
           ELSE
             IF(ISOM.EQ.2) THEN
              WRITE(IOUT,72) JPLAN,JFICH,DTHETAP,DPHIP,ECIN,SSETDIR_1,SS
     &ETDIF_1,SSETDIR_2,SSETDIF_2
              IF(ICHKDIR.EQ.2) THEN
               WRITE(IOUT,72) JPLAN,JFICH,DTHETAP,DPHIP,ECIN,SSET2DIR_1,
     &SSET2DIF_1,SSET2DIR_2,SSET2DIF_2
              ENDIF
             ELSE
              WRITE(IOUT,72) JPLAN,JEMET,DTHETAP,DPHIP,ECIN,SSETDIR_1,SS
     &ETDIF_1,SSETDIR_2,SSETDIF_2
              IF(ICHKDIR.EQ.2) THEN
               WRITE(IOUT,72) JPLAN,JEMET,DTHETAP,DPHIP,ECIN,SSET2DIR_1,
     &SSET2DIF_1,SSET2DIR_2,SSET2DIF_2
              ENDIF
             ENDIF
           ENDIF
 222       CONTINUE
C
C     End of the loop on the scanned angle
C
          ENDDO
C
   8      CONTINUE
C
C     End of the loop on the fixed angle
C
         ENDDO
C
C     End of the loop on the energy
C
         CLOSE(IUI6)
        ENDDO
C
   3    CONTINUE
C
C     End of the loop on the emitters
C
       ENDDO
C
       GO TO 1
   5   IPLAN=JPLAN-1
       IJK=IJK+1
       IF((IJK.EQ.1).AND.(IPRINT.GT.0)) THEN
         IF(I_TEST.NE.2) WRITE(IUO1,54) IPLAN
       ENDIF
   1   CONTINUE
C
C     End of the loop on the planes
C
      ENDDO
C
      IF(ABS(I_EXT).GE.1) CLOSE(IUI6)
      IF((ISOM.EQ.0).OR.(JFICH.EQ.NFICHLEC)) WRITE(IOUT,*)
      IF(SPECTRO.EQ.'APC') CLOSE(IOUT)
      IF(SPECTRO.EQ.'APC') GOTO 7
c      IF(((NEMET.GT.1).OR.(NPLAN.GT.1)).AND.(ISOM.EQ.0)) THEN
      IF(((NEMET.GT.1).OR.(NPLAN.GT.0)).AND.(ISOM.EQ.0)) THEN
        NP=0
        CALL TREAT_PHD(ISOM,NFICHLEC,JFICH,NP)
      ENDIF
      IF(I_EXT.EQ.2) THEN
         CALL WEIGHT_SUM(ISOM,I_EXT,0,1)
      ENDIF
      GOTO 7
   6  WRITE(IUO1,55)
C
   9  FORMAT(9(2X,I1),2X,I2)
  11  FORMAT(I4,2X,I4,2X,I4)
  12  FORMAT(2X,A3,11X,A13)
  13  FORMAT(6X,I1,1X,I3,2X,I4)
  14  FORMAT(6X,I1,1X,I3,3X,I3)
  19  FORMAT(2(2X,I1),1X,I2,6(2X,I1),2X,I2)
  20  FORMAT(2(5X,F6.2,2X,F6.2),2X,I1)
  21  FORMAT(10X,E12.6,3X,E12.6)
  22  FORMAT(16X,'INTERNAL CALCULATION OF MEAN SQUARE DISPLACEMENTS',/,2
     &5X,' BY DEBYE UNCORRELATED MODEL:',/)
  23  FORMAT(21X,'ATOM TYPE ',I5,'  MSD = ',F8.6,' ANG**2')
  51  FORMAT(/////,2X,'******* PLANE NUMBER ',I3,' DOES NOT CONTAIN ',
     *'ANY ABSORBER OF TYPE ',I2,' *******')
  52  FORMAT(/////,2X,'******* PLANE NUMBER ',I3,' POSITION OF ','THE AB
     &SORBER : (',F6.3,',',F6.3,',',F6.3,') *******',/,2X,'******* ',19X
     &,'THIS ABSORBER IS OF TYPE ',I2,20X,' *******')
  53  FORMAT(//,2X,'ORDER ',I2,'  TOTAL NUMBER OF PATHS     : ',F15.1,/,
     &10X,'  EFFECTIVE NUMBER OF PATHS : ',F15.1,/,10X,'  MINIMAL INTENS
     &ITY         : ',E12.6,2X,'No OF THE PATH : ',F15.1,               
     &          /,10X,'  MAXIMAL INTENSITY         : ',E12.6,2X,'No OF T
     &HE PATH : ',F15.1)
  54  FORMAT(//,7X,'DUE TO THE SIZE OF THE CLUSTER, THE SUMMATION',
     *' HAS BEEN TRUNCATED TO THE ',I2,' TH PLANE')
  55  FORMAT(///,12X,' <<<<<<<<<<  THIS VALUE OF ILPM IS NOT',
     *'AVAILABLE  >>>>>>>>>>')
  56  FORMAT(4X,'LATTICE PARAMETER A = ',F6.3,' ANGSTROEMS',4X,
     *'MEAN FREE PATH = ',F6.3,' * A',//)
  57  FORMAT(25X,'CLUSTER RADIUS = ',F6.3,' *A')
  58  FORMAT(//,2X,'ORDER ',I2,'  TOTAL NUMBER OF PATHS     : ',I10,/,10
     &X,'  EFFECTIVE NUMBER OF PATHS : ',I10,                /,10X,'  MI
     &NIMAL INTENSITY         : ',E12.6,2X,'No OF THE PATH : ',I10,     
     &                    /,10X,'  MAXIMAL INTENSITY         : ',E12.6, 
     &            2X,'No OF THE PATH : ',I10)
  59  FORMAT(//,15X,'THE SCATTERING DIRECTION IS GIVEN INSIDE ',
     *'THE CRYSTAL')
  60  FORMAT(7X,'THE POSITIONS OF THE ATOMS ARE GIVEN WITH RESPECT ',
     *'TO THE ABSORBER')
  61  FORMAT(///,4X,'..........  DIRECTION OF THE DETECTOR      : (',F6.
     &3,',',F6.3,',',F6.3,    ')  ..........',/,16X,'DIRECTION OF THE LI
     &GHT    ',     '     : (',F6.3,',',F6.3,',',F6.3,                  
     &                     ')',/,16X,'DIRECTION OF THE POLARIZATION  : (
     &',     F6.3,',',F6.3,',',F6.3,')',/,16X,'ANALYZER.LIGHT       ',' 
     &         :        ',F7.4)
  63  FORMAT(///,4X,'..........  DIRECTION OF THE DETECTOR      : (',F6.
     &3,',',F6.3,',',F6.3,    ')  ..........',/,16X,'DIRECTION OF THE LI
     &GHT    ',     '     : (',F6.3,',',F6.3,',',F6.3,')',/,16X,'ANALYZE
     &R.LIGHT               :        ',F7.4)
  65  FORMAT(////,3X,'++++++++++++++++++',9X,
     *'THETA = ',F6.2,' DEGREES',9X,'++++++++',
     *'++++++++++',///)
  66  FORMAT(////,3X,'++++++++++++++++++',9X,
     *'PHI = ',F6.2,' DEGREES',9X,'++++++++++',
     *'++++++++++',///)
  67  FORMAT(2X,I3,2X,I2,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
  68  FORMAT(10X,'  CUT-OFF INTENSITY       : ',E12.6)
  69  FORMAT(9X,I2,2X,E12.6,7X,E12.6,1X,F6.3,1X,10(I3,2X))
  70  FORMAT(2X,I2,2X,I10,7X,E12.6,2X,F6.3,7X,I2,7X,10(I3,2X))
  71  FORMAT(//,1X,'JDIF',4X,'No OF THE PATH',2X,'INTENSITY',3X,'LENGTH'
     &,4X,'ABSORBER',2X,'ORDER OF THE SCATTERERS',/)
  72  FORMAT(2X,I3,2X,I2,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6,2X,E1
     &2.6,2X,E12.6)
  74  FORMAT(10X,'<===== NUMBER OF PATHS TOO LARGE FOR PRINTING ','=====
     &>')
  76  FORMAT(2X,I2,2X,E12.6,7X,E12.6,2X,F6.3,7X,I2,7X,10(I3,2X))
  77  FORMAT('   ')
  79  FORMAT(2X,I3,2X,I2,2X,I4,2X,I4,2X,I4)
  80  FORMAT(///)
  81  FORMAT(//,1X,'RANK',1X,'ORDER',4X,'No PATH',3X,'INTENSITY',3X,'LEN
     &GTH',4X,'ABS',3X,'ORDER OF THE SCATTERERS',/)
  82  FORMAT(I3,4X,I2,1X,E12.6,3X,E12.6,2X,F6.3,4X,I2,4X,10(I3,1X))
  83  FORMAT(I3,4X,I2,1X,I10,3X,E12.6,2X,F6.3,4X,I2,4X,10(I3,1X))
  84  FORMAT(/////,18X,'THE ',I3,' MORE INTENSE PATHS BY DECREASING',' O
     &RDER :',/,24X,'(THE LENGTH IS GIVEN IN UNITS ','OF A)')
  85  FORMAT(/////,25X,' PATHS USED IN THE CALCULATION : ',/,24X,'(THE L
     &ENGTH IS GIVEN IN UNITS OF A)')
  86  FORMAT(2X,I3,1X,I4,5X,F8.3,3X,F8.3,3X,E12.6)
  87  FORMAT(2X,I2,2X,I3,2X,I2,2X,I3,2X,I3,2X,I3,2X,I1,2X,I2,2X,I2,2X,E1
     &2.6,2X,E12.6,2X,E12.6,2X,E12.6)
  88  FORMAT(/,19X,'TILTED THETA =',F6.2,5X,'TILTED PHI =', F6.2)
  89  FORMAT(/,4X,'..........................................','........
     &.............................')
C
   7  CALL PATH_SAMPLING_SUMMARY()
      RETURN
C
      END
