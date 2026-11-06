C
C
C           ************************************************************
C           * ******************************************************** *
C           * *                                                      * *
C           * *         MULTIPLE-SCATTERING SPIN-INDEPENDENT         * *
C           * *            PHOTOELECTRON DIFFRACTION CODE            * *
C           * *            BASED ON CORRELATION EXPANSION            * *
C           * *                                                      * *
C           * ******************************************************** *
C           ************************************************************
C
C
C
C
C  Written by D. Sebilleau, Groupe Theorie,
C                           Departement Materiaux-Nanosciences,
C                           Institut de Physique de Rennes,
C                           UMR CNRS-Universite 6251,
C                           Universite de Rennes-1,
C                           35042 Rennes-Cedex,
C                           France
C
C  Contributions : M. Gavaza, H.-F. Zhao, K. Hatada
C
C-----------------------------------------------------------------------
C
C     As a general rule in this code, although there might be a few
C     exceptions (...), a variable whose name starts with a 'I' is a 
C     switch, with a 'J' is a loop index and with a 'N' is a number.
C
C     The main subroutines are :
C
C                * PHDDIF    : computes the photoelectron diffraction 
C                              formula
C
C                * LEDDIF    : computes the low-energy electron  
C                              diffraction formula
C
C                * XASDIF    : computes the EXAFS or XANES formula
C                              depending on the energy 
C
C                * AEDDIF    : computes the Auger electron diffraction
C                              formula 
C
C
C     A subroutine called NAME_A is the Auger equivalent of subroutine 
C     NAME. The essentail difference between NAME and NAME_A is that 
C     they do not contain the same arrays.
C
C     Always remember, when changing the input data file, to keep the 
C     format. The rule here is that the last digit of any integer or 
C     character data must correspond to the tab (+) while for real data,     
C     the tab precedes the point. 
C
C     Do not forget, before submitting a calculation, to check the 
C     consistency of the input data with the corresponding maximal
C     values in the include file. 
C
C-----------------------------------------------------------------------
C
C     Please report any bug or problem to me at : 
C
C                      didier.sebilleau@univ-rennes1.fr
C
C
C
C  Last modified : 10 Jan 2016
C
C=======================================================================
C
      SUBROUTINE MAIN_PHD_NS_CE()
C
C  This routine reads the various input files and calls the subroutine
C               performing the requested calculation
C
      USE DIM_MOD
      USE ADSORB_MOD
      USE APPROX_MOD
      USE ATOMS_MOD
      USE AUGER_MOD
      USE BASES_MOD
      USE CLUSLIM_MOD
      USE COOR_MOD
      USE DEBWAL_MOD
      USE INDAT_MOD
      USE INIT_A_MOD
      USE INIT_L_MOD
      USE INIT_J_MOD
      USE INIT_M_MOD
      USE INFILES_MOD
      USE INUNITS_MOD
      USE LIMAMA_MOD
      USE LPMOY_MOD
      USE MASSAT_MOD
      USE MILLER_MOD
      USE OUTUNITS_MOD
      USE PARCAL_MOD
      USE PARCAL_A_MOD
      USE RELADS_MOD
      USE RELAX_MOD
      USE RESEAU_MOD
      USE SPIN_MOD
      USE TESTS_MOD
      USE TRANS_MOD
      USE TL_AED_MOD
      USE TYPCAL_MOD
      USE TYPCAL_A_MOD
      USE TYPEM_MOD
      USE TYPEXP_MOD
      USE VALIN_MOD
      USE XMRHO_MOD
C
      DIMENSION VEC(3,3),VB1(3),VB2(3),VB3(3),VBS(3)
      DIMENSION ROT(3,3),EMET(3)
      DIMENSION VAL2(NATCLU_M)
      DIMENSION IRE(NATCLU_M,2)
      DIMENSION REL(NATCLU_M),RHOT(NATM)
      DIMENSION ATOME(3,NATCLU_M),COORD(3,NATCLU_M)
      DIMENSION NTYP(NATCLU_M),NATYP_OLD(NATM)
      DIMENSION LMAX_TMP(NATM,NE_M),DIST12(NATCLU_M,NATCLU_M)
      DIMENSION IBWD_TMP(NATP_M),RTHFWD_TMP(NATP_M),RTHBWD_TMP(NATP_M)
      DIMENSION UJ2_TMP(NATM),RHOT_TMP(NATM),XMT_TMP(NATM)
C
      COMPLEX TLSTAR
      COMPLEX RHOR(NE_M,NATM,0:18,5,NSPIN2_M)
      COMPLEX TLSTAR_A
      COMPLEX RHOR_A(0:NT_M,NATM,0:40,2,NSPIN2_M),RAD_D,RAD_E
      COMPLEX RHOR1STAR,RHOR2STAR,RHOR3STAR,RHOR4STAR,RHOR5STAR
C
      INTEGER INV(2)
C
      CHARACTER RIEN
      CHARACTER*1 B
      CHARACTER*2 R
C
C
C
C
C
C
      CHARACTER*30 TUNIT,DUMMY
C
      DATA PI,BOHR,SMALL/3.141593,0.529177,0.001/
      DATA INV /0,0/
C
      LE_MAX=0
C
C!      READ(*,776) NFICHLEC
C!      READ(*,776) ICOM
C!      DO JF=1,NFICHLEC
C!        READ(*,777) INDATA(JF)
C!      ENDDO
C
C..........  Loop on the data files  ..........
C
      NFICHLEC=1
      ICOM = 5
      DO JFICH=1,NFICHLEC
C!      OPEN(UNIT=ICOM, FILE=INDATA(JFICH), STATUS='OLD')
      OPEN(UNIT=ICOM, FILE='../input/spec.dat', STATUS='OLD')
      CALL READ_DATA(ICOM,NFICHLEC,JFICH,ITRTL,*2,*1,*55,*74,*99,*504,*5
     &20,*540,*550,*570,*580,*590,*630)
C
C..........  Atomic case index  ..........
C
      I_AT=0
      IF((SPECTRO.EQ.'PHD').AND.(I_TEST.EQ.2)) I_AT=1
      IF((SPECTRO.EQ.'LED').AND.(I_TEST.EQ.2)) I_AT=1
      IF((SPECTRO.EQ.'AED').AND.(I_TEST_A.EQ.2)) I_AT=1
      IF((SPECTRO.EQ.'XAS').AND.(I_TEST.EQ.2)) I_AT=1
      IF(SPECTRO.EQ.'APC') THEN 
        IF((I_TEST.EQ.2).AND.(I_TEST_A.EQ.2)) I_AT=1
      ENDIF
C
      IF(IBAS.EQ.1) THEN
        IF(ITEST.EQ.0) THEN
          NEQ=(2*NIV+1)**3
        ELSE
          NEQ=(2*NIV+3)**3
        ENDIF
        IF(NEQ*NATP_M.GT.NATCLU_M) GOTO 518
      ENDIF
C
      IF(SPECTRO.EQ.'APC') THEN
        N_EL=2
      ELSE
        N_EL=1
      ENDIF
      IF((INTERACT.EQ.'COULOMB').OR.(INTERACT.EQ.'DIPCOUL')) THEN
        IF(I_MULT.EQ.0) THEN
          LE_MIN=ABS(LI_C-ABS(LI_I-LI_A))
          LE_MAX=LI_C+LI_A+LI_I
        ELSE
          LE_MIN=ABS(LI_C-L_MUL)
          LE_MAX=LI_C+L_MUL
        ENDIF
      ENDIF
C
C..........  Test of the dimensions against the input values  ..........
C
      IF(NO.GT.NO_ST_M) GOTO 600
      IF(LE_MAX.GT.LI_M) GOTO 620
C
      OPEN(UNIT=IUI2, FILE=INFILE2, STATUS='OLD')
      OPEN(UNIT=IUI3, FILE=INFILE3, STATUS='OLD')
      IF(INTERACT.EQ.'DIPCOUL') THEN
        OPEN(UNIT=IUI7, FILE=INFILE7, STATUS='OLD')
        OPEN(UNIT=IUI8, FILE=INFILE8, STATUS='OLD')
      ENDIF
C
C..........  Reading of the TL and radial matrix elements files  ..........
C..........      (dipolar excitation or no excitation case)       ..........
C
      IF(INTERACT.NE.'COULOMB') THEN
        IF(SPECTRO.EQ.'APC') WRITE(IUO1,418)
        READ(IUI2,3) NAT1,NE1,ITL,IPOTC,LMAX_MODE
        IF(ISPIN.EQ.0) THEN
          IF(NAT1.EQ.1) THEN
            WRITE(IUO1,561)
          ELSE
            WRITE(IUO1,560) NAT1
          ENDIF
        ENDIF
        IF((ITL.EQ.1).AND.(ISPIN.EQ.1)) THEN
          READ(IUI2,530) E_MIN,E_MAX,DE
        ENDIF
        IF((ISPIN.EQ.0).AND.(ITL.EQ.0)) THEN
          NLG=INT(NAT1-0.0001)/4 +1
          DO NN=1,NLG
            NRL=4*NN
            JD=4*(NN-1)+1
            IF(NN.EQ.NLG) NRL=NAT1
            READ(IUI2,555) (LMAX(JAT,1),JAT=JD,NRL)
            WRITE(IUO1,556) (LMAX(JAT,1),JAT=JD,NRL)
          ENDDO
C
C   Temporary storage of LMAX. Waiting for a version of PHAGEN
C           with LMAX dependent on the energy
C
          DO JE=1,NE
            DO JAT=1,NAT1
              LMAX(JAT,JE)=LMAX(JAT,1)
            ENDDO 
          ENDDO
C            
          NL1=1
          DO JAT=1,NAT1
            NL1=MAX0(NL1,LMAX(JAT,1)+1)
          ENDDO
          IF(NL1.GT.NL_M) GOTO 184
        ENDIF
        IF(ITL.EQ.0) READ(IUI3,101) NATR,NER
        IF(ISPIN.EQ.1) THEN
          READ(IUI3,106) L_IN,NATR,NER
          IF(LI.NE.L_IN) GOTO 606
        ENDIF
        NAT2=NAT+NATA
        IF((NAT1.NE.NAT2).OR.(NE1.NE.NE)) GOTO 180
        IF((ITL.EQ.0).AND.((NATR.NE.NAT2).OR.(NER.NE.NE))) GOTO 182
C
C..........  DL generated by MUFPOT and RHOR given  ..........
C..........     by S. M. Goldberg, C. S. Fadley     ..........
C..........     and S. Kono, J. Electron Spectr.    ..........
C..........      Relat. Phenom. 21, 285 (1981)      ..........
C
        IF(ITL.EQ.0) THEN
          DO JAT=1,NAT2     
            IF((INITL.NE.0).AND.(IFTHET.NE.1)) THEN
              READ(IUI3,102) RIEN  
              READ(IUI3,102) RIEN   
              READ(IUI3,102) RIEN   
            ENDIF
            DO JE=1,NE
              IF((IFTHET.EQ.1).OR.(INITL.EQ.0)) GOTO 121
              READ(IUI3,103) ENERGIE
              READ(IUI3,102) RIEN
              READ(IUI3,102) RIEN 
              READ(IUI3,102) RIEN  
 121          CONTINUE
              DO L=0,LMAX(JAT,JE)
                READ(IUI2,7) VK(JE),TL(L,1,JAT,JE)
                TL(L,1,JAT,JE)=CSIN(TL(L,1,JAT,JE))*CEXP((0.,1.)*
     1                         TL(L,1,JAT,JE))
              ENDDO
              IF((IFTHET.EQ.1).OR.(INITL.EQ.0)) GOTO 5
              DO LL=1,18
                READ(IUI3,104) RH1,RH2,DEF1,DEF2
                RHOR(JE,JAT,LL,1,1)=CMPLX(RH1)
                RHOR(JE,JAT,LL,2,1)=CMPLX(RH2)
                DLT(JE,JAT,LL,1)=CMPLX(DEF1)
                DLT(JE,JAT,LL,2)=CMPLX(DEF2)
              ENDDO
  5           CONTINUE          
            ENDDO            
          ENDDO              
        ELSE
C
C..........  TL and RHOR calculated by PHAGEN  ..........
C
          DO JE=1,NE
            NLG=INT(NAT2-0.0001)/4 +1
            IF(NE.GT.1) WRITE(IUO1,563) JE
            DO NN=1,NLG
              NRL=4*NN
              JD=4*(NN-1)+1
              IF(NN.EQ.NLG) NRL=NAT2
              READ(IUI2,555) (LMAX(JAT,JE),JAT=JD,NRL)
              WRITE(IUO1,556) (LMAX(JAT,JE),JAT=JD,NRL)
            ENDDO
            NL1=1
            DO JAT=1,NAT2
              NL1=MAX0(NL1,LMAX(JAT,1)+1)
            ENDDO
            IF(NL1.GT.NL_M) GOTO 184
            DO JAT=1,NAT2
              READ(IUI2,*) DUMMY
              DO L=0,LMAX(JAT,JE)
                IF(LMAX_MODE.EQ.0) THEN
                  READ(IUI2,9) VK(JE),TLSTAR
                ELSE
                  READ(IUI2,9) VK(JE),TLSTAR
                ENDIF
                TL(L,1,JAT,JE)=CONJG(TLSTAR)
                VK(JE)=CONJG(VK(JE))
              ENDDO
            ENDDO
C
            IF((IFTHET.EQ.1).OR.(INITL.EQ.0)) GOTO 333
            IF(JE.EQ.1) THEN
              DO JDUM=1,7
                 READ(IUI3,102) RIEN
              ENDDO
            ENDIF
C
C  Reading or regular (RHOR) and irregular (RHOI) radial integrals
C
C                  1-2 : dipole terms
C                  3-5 : quadrupole terms
C
            DO JEMET=1,NEMET
C
              JM=IEMET(JEMET)
              READ(IUI3,105) RHOR1STAR,RHOR2STAR,RHOR3STAR,RHOR4STAR,
     1                       RHOR5STAR
              RHOR(JE,JM,NNL,1,1)=CONJG(RHOR1STAR)
              RHOR(JE,JM,NNL,2,1)=CONJG(RHOR2STAR)    
              RHOR(JE,JM,NNL,3,1)=CONJG(RHOR3STAR)    
              RHOR(JE,JM,NNL,4,1)=CONJG(RHOR4STAR)    
              RHOR(JE,JM,NNL,5,1)=CONJG(RHOR5STAR)
C    
            ENDDO
C
  333       VK(JE)=VK(JE)*A
            VK2(JE)=CABS(VK(JE)*VK(JE))
          ENDDO
        ENDIF 
C
        CLOSE(IUI2)
        CLOSE(IUI3)
C 
C.......... Suppression of possible zeros in the TL array  ..........
C..........  (in case of the use of matrix inversion and   ..........
C..........             for energy variations)             ..........
C
        IF((ISPIN.EQ.0).AND.(ITL.EQ.1).AND.(LMAX_MODE.NE.0)) THEN
           CALL SUP_ZEROS(TL,LMAX,NE,NAT2,IUO1,ITRTL)
        ENDIF

      ENDIF
C
C..........  Reading of the TL and radial matrix elements files  ..........
C..........             (Coulomb excitation case)                ..........
C
      IF((INTERACT.EQ.'COULOMB').OR.(INTERACT.EQ.'DIPCOUL')) THEN
        IERR=0
        IF(INTERACT.EQ.'COULOMB') THEN
          IRD1=IUI2
          IRD2=IUI3
        ELSEIF(INTERACT.EQ.'DIPCOUL') THEN
          IRD1=IUI7
          IRD2=IUI8
        ENDIF
        IF(SPECTRO.EQ.'APC') WRITE(IUO1,419)
        READ(IRD1,3) NAT1_A,NE1_A,ITL_A,IPOTC_A,LMAX_MODE_A
        IF(ISPIN.EQ.0) THEN
          IF(NAT1_A.EQ.1) THEN
            WRITE(IUO1,561)
          ELSE
            WRITE(IUO1,560) NAT1_A
          ENDIF
        ENDIF
        IF((ITL_A.EQ.1).AND.(ISPIN.EQ.1)) THEN
          READ(IRD1,530) E_MIN_A,E_MAX_A,DE_A
        ENDIF
        IF(ITL_A.EQ.1) THEN
          READ(IRD2,107) LI_C2,LI_I2,LI_A2
          READ(IRD2,117) LE_MIN1,N_CHANNEL
          LE_MAX1=LE_MIN1+N_CHANNEL-1
          IF(I_TEST_A.NE.1) THEN
            IF((LE_MIN.NE.LE_MIN1).OR.(LE_MAX.NE.LE_MAX1)) GOTO 610
          ELSE
            LI_C2=0
            LI_I2=1
            LI_A2=0
            LE_MIN1=1
            N_CHANNEL=1
          ENDIF
        ENDIF
        IF((ISPIN.EQ.0).AND.(ITL_A.EQ.0)) THEN
          NLG=INT(NAT1_A-0.0001)/4 +1
          DO NN=1,NLG
            NRL=4*NN
            JD=4*(NN-1)+1
            IF(NN.EQ.NLG) NRL=NAT1_A
            READ(IRD1,555) (LMAX_A(JAT,1),JAT=JD,NRL)
            WRITE(IUO1,556) (LMAX_A(JAT,1),JAT=JD,NRL)
          ENDDO
C
C   Temporary storage of LMAX_A. Waiting for a version of PHAGEN
C           with LMAX_A dependent on the energy
C
          DO JE=1,NE1_A
            DO JAT=1,NAT1_A
              LMAX_A(JAT,JE)=LMAX_A(JAT,1)
            ENDDO 
          ENDDO
C            
          NL1_A=1
          DO JAT=1,NAT1_A
            NL1_A=MAX0(NL1_A,LMAX_A(JAT,1)+1)
          ENDDO
          IF(NL1_A.GT.NL_M) GOTO 184
        ENDIF
        IF(ITL_A.EQ.0) READ(IRD2,101) NATR_A,NER_A
        IF(ISPIN.EQ.1) THEN
          READ(IRD2,106) L_IN_A,NATR_A,NER_A
          IF(LI_C.NE.L_IN_A) GOTO 606
        ENDIF
        NAT2_A=NAT+NATA
        NAT2=NAT2_A
        IF((NAT1_A.NE.NAT2_A).OR.(NE1_A.NE.NE_A)) GOTO 180
        IF((ITL_A.EQ.0).AND.((NATR_A.NE.NAT2_A).OR.(NER_A.NE.NE))) 
     1    GOTO 182
C
C..........  DL generated by MUFPOT and RHOR given  ..........
C..........     by S. M. Goldberg, C. S. Fadley     ..........
C..........     and S. Kono, J. Electron Spectr.    ..........
C..........      Relat. Phenom. 21, 285 (1981)      ..........
C
        IF(ITL_A.EQ.0) THEN
          CONTINUE
        ELSE
C
C..........  TL_A and RHOR_A calculated by PHAGEN  ..........
C
          DO JE=1,NE_A
            NLG=INT(NAT2_A-0.0001)/4 +1
            IF(NE_A.GT.1) WRITE(IUO1,563) JE
            DO NN=1,NLG
              NRL=4*NN
              JD=4*(NN-1)+1
              IF(NN.EQ.NLG) NRL=NAT2_A
              READ(IRD1,555) (LMAX_A(JAT,JE),JAT=JD,NRL)
              WRITE(IUO1,556) (LMAX_A(JAT,JE),JAT=JD,NRL)
            ENDDO
            DO JAT=1,NAT2_A
              READ(IRD1,*) DUMMY
              DO L=0,LMAX_A(JAT,JE)
                IF(LMAX_MODE_A.EQ.0) THEN
                  READ(IRD1,9) VK_A(JE),TLSTAR
                ELSE
                  READ(IRD1,7) VK_A(JE),TLSTAR
                ENDIF
                TL_A(L,1,JAT,JE)=CONJG(TLSTAR)
                VK_A(JE)=CONJG(VK_A(JE))
              ENDDO
            ENDDO
C
            IF(IFTHET_A.EQ.1) GOTO 331
            DO LE=LE_MIN,LE_MAX
              DO JEMET=1,NEMET
                JM=IEMET(JEMET)
                READ(IRD2,109) L_E,LB_MIN,LB_MAX
                IF(I_TEST_A.EQ.1) THEN
                  L_E=1
                  LB_MIN=0
                  LB_MAX=1
                ENDIF
                IF(LE.NE.L_E) IERR=1
                L_BOUNDS(L_E,1)=LB_MIN
                L_BOUNDS(L_E,2)=LB_MAX
                DO LB=LB_MIN,LB_MAX
                  READ(IRD2,108) L_A,RAD_D,RAD_E
                  RHOR_A(LE,JM,L_A,1,1)=RAD_D
                  RHOR_A(LE,JM,L_A,2,1)=RAD_E
                  IF(I_TEST_A.EQ.1) THEN
                    IF(LB.EQ.LB_MIN) THEN
                      RHOR_A(LE,JM,L_A,1,1)=(0.0,0.0)
                      RHOR_A(LE,JM,L_A,2,1)=(1.0,0.0)
                    ELSEIF(LB.EQ.LB_MAX) THEN
                      RHOR_A(LE,JM,L_A,1,1)=(1.0,0.0)
                      RHOR_A(LE,JM,L_A,2,1)=(0.0,0.0)
                    ENDIF
                  ENDIF
                ENDDO 
              ENDDO
            ENDDO
  331       VK_A(JE)=VK_A(JE)*A
            VK2_A(JE)=CABS(VK_A(JE)*VK_A(JE))
          ENDDO
        ENDIF 
C
        CLOSE(IRD1)
        CLOSE(IRD2)
C 
C.......... Suppression of possible zeros in the TL array  ..........
C..........  (in case of the use of matrix inversion and   ..........
C..........             for energy variations)             ..........
C
        IF((ISPIN.EQ.0).AND.(ITL_A.EQ.1).AND.(LMAX_MODE_A.NE.0)) THEN
           CALL SUP_ZEROS(TL_A,LMAX_A,NE_A,NAT2_A,IUO1,ITRTL)
        ENDIF
        IF(SPECTRO.EQ.'APC') WRITE(IUO1,420)
C
      ENDIF
C
C..........    Checking maximum value for l_max     ..........
C..........    and storage of Gaunt coefficients    ..........
C
      LM_PE=0
      DO JAT=1,NAT2
        DO JE=1,NE
          LM_PE=MAX(LM_PE,LMAX(JAT,JE))
        ENDDO
      ENDDO
C
      LM_AE=0
      DO JAT=1,NAT2_A
        DO JE=1,NE_A
          LM_AE=MAX(LM_AE,LMAX_A(JAT,JE))
        ENDDO
      ENDDO
C
      LM_PA=MAX(LM_PE,LM_AE)
      CALL GAUNT_ST(LM_PA)
      CALL COEFPQ(MAX(NAT2,NAT2_A),NDIF)
C
C..........  Check of the consistency of the two TL and radial ..........
C..........           matrix elements for APECS                ..........
C
      IF(SPECTRO.EQ.'APC') THEN
C
        I_TL_FILE=0
        I_RD_FILE=0
C
        IF(NAT1.NE.NAT1_A) I_TL_FILE=1
        IF(NE1.NE.NE1_A) I_TL_FILE=1
        IF(ITL.NE.ITL_A) I_TL_FILE=1
        IF(IPOTC.NE.IPOTC_A) I_TL_FILE=1
C
        IF(LI_C.NE.LI_C2) I_RD_FILE=1
        IF(LI_I.NE.LI_I2) I_RD_FILE=1
        IF(LI_A.NE.LI_A2) I_RD_FILE=1
C
        IF(I_TL_FILE.EQ.1) GOTO 608
        IF(I_RD_FILE.EQ.1) GOTO 610
        IF(IERR.EQ.1) GOTO 610
C
      ENDIF
C
C..........  Calculation of the scattering factor (only)  ..........
C
      IF((IFTHET.EQ.0).AND.(IFTHET_A.EQ.0)) GO TO 8   
      IF(IFTHET.EQ.1) THEN 
        CALL PLOTFD(A,LMAX,ITL,NL1,NAT2,NE)
      ELSEIF(IFTHET_A.EQ.1) THEN
c        CALL PLOTFD_A(A,LMAX_A,ITL_A,NL1_A,NAT2_A,NE_A)
      ENDIF
      WRITE(IUO1,57)
CST   STOP
      GO TO 999
C
   8  IF(IBAS.EQ.0) THEN
C
C...............  Reading of an external cluster  ...............
C
C
C  Cluster originating from CLUSTER_NEW.F : IPHA=0
C  Cluster originating from PHAGEN_NEW.F  : IPHA=1 (atomic units), IPHA=2 (angstroems)
C  Other cluster                          : the first line must be text; then
C                                            free format : Atomic number,X,Y,Z,number
C                                            of the corresponding prototypical atom ;
C                                            All atoms corresponding to the same 
C                                            prototypical atom must follow each other.
C                                            Moreover, the blocks of equivalent atoms
C                                            must be ordered by increasing number of 
C                                            prototypical atom.
C
        VALZ_MIN=1000.0
        VALZ_MAX=-1000.0
C
        OPEN(UNIT=IUI4, FILE=INFILE4, STATUS='OLD')
        READ(IUI4,778,ERR=892) IPHA
        GOTO 893
  892   IPHA=3
        IF(UNIT.EQ.'ANG') THEN
          CUNIT=1./A
          TUNIT='ANGSTROEMS'
        ELSEIF(UNIT.EQ.'LPU') THEN
          CUNIT=1.
          TUNIT='UNITS OF THE LATTICE PARAMETER'
        ELSEIF(UNIT.EQ.'ATU') THEN
          CUNIT=BOHR/A
          TUNIT='ATOMIC UNITS'
        ELSE
          GOTO 890
        ENDIF
  893   NATCLU=0
        DO JAT=1,NAT2
          NATYP(JAT)=0
        ENDDO
        IF(IPHA.EQ.0) THEN
          CUNIT=1.
          TUNIT='UNITS OF THE LATTICE PARAMETER'
        ELSEIF(IPHA.EQ.1) THEN
          CUNIT=BOHR/A
          TUNIT='ATOMIC UNITS'
          IEMET(1)=1
        ELSEIF(IPHA.EQ.2) THEN
          CUNIT=1./A
          TUNIT='ANGSTROEMS'
          IEMET(1)=1
        ENDIF
        IF(IPRINT.EQ.2) THEN
          IF(I_AT.NE.1) THEN
            WRITE(IUO1,558) IUI4,TUNIT
            IF(IPHA.EQ.3) WRITE(IUO1,549) 
          ENDIF
        ENDIF
        JATM=0
        DO JLINE=1,10000
          IF(IPHA.EQ.0) THEN
            READ(IUI4,125,END=780) R,NN,X,Y,Z,JAT
          ELSEIF(IPHA.EQ.1) THEN
            READ(IUI4,779,END=780) R,NN,X,Y,Z,JAT
          ELSEIF(IPHA.EQ.2) THEN
            READ(IUI4,779,END=780) R,NN,X,Y,Z,JAT
          ELSEIF(IPHA.EQ.3) THEN
            READ(IUI4,*,END=780) NN,X,Y,Z,JAT
          ENDIF
          JATM=MAX0(JAT,JATM)
          NATCLU=NATCLU+1
          IF(IPHA.NE.3) THEN
            CHEM(JAT)=R
          ELSE
            CHEM(JAT)='XX'
          ENDIF
          NZAT(JAT)=NN
          NATYP(JAT)=NATYP(JAT)+1
          COORD(1,NATCLU)=X*CUNIT
          COORD(2,NATCLU)=Y*CUNIT
          COORD(3,NATCLU)=Z*CUNIT
          VALZ(NATCLU)=Z*CUNIT
          IF((IPRINT.GE.2).AND.(I_AT.EQ.0)) THEN
            WRITE(IUO1,557) NATCLU,COORD(1,NATCLU),COORD(2,NATCLU),
     1                      COORD(3,NATCLU),JAT,NATYP(JAT),CHEM(JAT)
          ENDIF
        ENDDO
 780    NBZ=NATCLU
        IF(JATM.NE.NAT) GOTO 514
        CLOSE(IUI4)
C
        IF(NATCLU.GT.NATCLU_M) GOTO 510
        DO JA1=1,NATCLU
          DO JA2=1,NATCLU
            DIST12(JA1,JA2)=SQRT((COORD(1,JA1)-COORD(1,JA2))**2
     1                          +(COORD(2,JA1)-COORD(2,JA2))**2
     2                          +(COORD(3,JA1)-COORD(3,JA2))**2)
            IF((JA2.GT.JA1).AND.(DIST12(JA1,JA2).LT.0.001)) GOTO 895
          ENDDO
        ENDDO
C        
        D_UP=VALZ_MAX-VALZ(1)
        D_DO=VALZ(1)-VALZ_MIN
        IF((D_DO.LE.D_UP).AND.(I_GR.EQ.2)) THEN
          I_INV=1
        ELSE
          I_INV=0
        ENDIF
      ELSE
C
C...............  Construction of an internal cluster  ...............
C
        CALL BASE
        CALL ROTBAS(ROT)                                                  
        IF(IVG0.EQ.2) THEN
          NMAX=NIV+1
        ELSE
          NMAX=(2*NIV+1)**3
        ENDIF
        IF((IPRINT.EQ.2).AND.(IVG0.LE.1)) THEN 
          WRITE(IUO1,37)                                                        
          WRITE(IUO1,38) NIV 
          DO NUM=1,NMAX
            CALL NUMAT(NUM,NIV,IA,IB,IC)                              
            WRITE(IUO1,17) NUM,IA,IB,IC                                        
          ENDDO                                                        
          WRITE(IUO1,39)                                                        
        ENDIF    
        CALL AMAS(NIV,ATOME,COORD,VALZ,IESURF,COUPUR,ROT,
     1            IRE,NATYP,NBZ,NAT2,NCOUCH,NMAX)   
        IF((IREL.GE.1).OR.(NRELA.GT.0)) THEN                              
          CALL RELA(NBZ,NPLAN,NAT2,VALZ,VAL2,VAL,COORD,NATYP,REL,
     1              NCOUCH)  
          IF(IREL.EQ.1) THEN
            DO JP=1,NPLAN
              VAL(JP)=VAL2(JP)
            ENDDO
          ENDIF  
        ENDIF                                 
      ENDIF
C 
C  Storage of the extremal values of x and y for each plane. They define
C    the exterior of the cluster when a new cluster has to be build to 
C    support a point-group
C
      IF(I_GR.GE.1) THEN     
        IF((IREL.EQ.0).OR.(IBAS.EQ.0)) THEN
          CALL ORDRE(NBZ,VALZ,NPLAN,VAL)                                  
          WRITE(IUO1,50) NPLAN                                               
          DO K=1,NPLAN 
            WRITE(IUO1,29) K,VAL(K)
            X_MAX(K)=0.                                            
            X_MIN(K)=0.                                            
            Y_MAX(K)=0.                                            
            Y_MIN(K)=0.                                            
          ENDDO                                                        
        ENDIF
        DO JAT=1,NATCLU
          X=COORD(1,JAT)
          Y=COORD(2,JAT)
          Z=COORD(3,JAT)
          DO JPLAN=1,NPLAN
            IF(ABS(Z-VAL(JPLAN)).LT.SMALL) THEN
              X_MAX(JPLAN)=MAX(X,X_MAX(JPLAN))
              X_MIN(JPLAN)=MIN(X,X_MIN(JPLAN))
              Y_MAX(JPLAN)=MAX(Y,Y_MAX(JPLAN))
              Y_MIN(JPLAN)=MIN(Y,Y_MIN(JPLAN))
            ENDIF
          ENDDO  
        ENDDO
      ENDIF
C
C  Instead of the symmetrization of the cluster (this version only)
C 
      N_PROT=NAT
      NAT_ST=0
      DO JTYP=1,JATM
       NB_AT=NATYP(JTYP)
       IF(NB_AT.GT.NAT_EQ_M) GOTO 614
       DO JA=1,NB_AT
         NAT_ST=NAT_ST+1
         NCORR(JA,JTYP)=NAT_ST
       ENDDO
      ENDDO
      DO JC=1,3
        DO JA=1,NATCLU
          SYM_AT(JC,JA)=COORD(JC,JA)
        ENDDO
      ENDDO
C
C  Checking surface-like atoms for mean square displacements
C                    calculations
C
      CALL CHECK_VIB(NAT2)
C
C..........  Set up of the variables used for an internal  ..........
C..........   calculation of the mean free path and/or of  ..........
C..........         the mean square displacements          ..........
C
      IF((IDCM.EQ.1).OR.(ILPM.EQ.1)) THEN
        DO JTYP=1,NAT2
          XMT(JTYP)=XMAT(NZAT(JTYP))
          RHOT(JTYP)=RHOAT(NZAT(JTYP))
        ENDDO
        XMTA=XMT(1)
        RHOTA=RHOT(1)
        NZA=NZAT(1)
      ENDIF
      IF(IDCM.GT.0) THEN   
        CALL CHNOT(3,VECBAS,VEC)
        DO J=1,3
          VB1(J)=VEC(J,1)
          VB2(J)=VEC(J,2)
          VB3(J)=VEC(J,3)
        ENDDO
        CPR=1.
        CALL PRVECT(VB2,VB3,VBS,CPR)
        VM=PRSCAL(VB1,VBS)
        QD=(6.*PI*PI*NAT/VM)**(1./3.)
      ENDIF
C
C..........  Writing of the contents of the cluster,  ..........
C..........  of the position of the different planes  ..........
C..........   and of their respective absorbers in    ..........
C..........          the control file IUO1            ..........
C
      IF(I_AT.EQ.1) GOTO 153
      IF((IPRINT.EQ.2).AND.(IBAS.GT.0)) THEN  
        WRITE(IUO1,40)
        NCA=0                                                        
        DO J=1,NAT                                                   
          DO I=1,NMAX
            NCA=NCA+1
            WRITE(IUO1,20) J,I                                                 
            WRITE(IUO1,21) (ATOME(L,NCA),L=1,3)                                
            K=IRE(NCA,1)                                                      
            IF(K.EQ.0) THEN                                                 
              WRITE(IUO1,22)                                                      
            ELSE                                                            
              WRITE(IUO1,23) (COORD(L,K),L=1,3),IRE(NCA,2)                     
            ENDIF                                                           
          ENDDO                                                        
        ENDDO                                                        
        WRITE(IUO1,41)                                                        
      ENDIF 
      IF(IBAS.EQ.1) THEN                                                            
        WRITE(IUO1,24) 
        NATCLU=0                                                         
        DO I=1,NAT                                                     
          NN=NATYP(I) 
          NATCLU=NATCLU+NATYP(I)                                                      
          WRITE(IUO1,26) NN,I                                                  
        ENDDO
        IF(IADS.EQ.1) NATCLU=NATCLU+NADS1+NADS2+NADS3
        WRITE(IUO1,782) NATCLU
        IF(NATCLU.GT.NATCLU_M) GOTO 516
        IF(IPRINT.EQ.3) WRITE(IUO1,559)
        IF(IPRINT.EQ.3) THEN
          NBTA=0
          DO JT=1,NAT2
            NBJT=NATYP(JT)
            DO JN=1,NBJT
              NBTA=NBTA+1
              WRITE(IUO1,557) NBTA,COORD(1,NBTA),COORD(2,NBTA),
     1                        COORD(3,NBTA),JT,JN,CHEM(JT)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
 153  IF((ITEST.EQ.1).AND.(IBAS.GT.0)) THEN
        CALL TEST(NIV,ROT,NATYP,NBZ,NAT2,IESURF,COUPUR,*56) 
      ENDIF
      IF((IREL.EQ.0).OR.(IBAS.EQ.0)) THEN
        CALL ORDRE(NBZ,VALZ,NPLAN,VAL)                                  
        IF(I_AT.EQ.0) WRITE(IUO1,50) NPLAN                                               
        DO K=1,NPLAN 
          IF(I_AT.EQ.0) WRITE(IUO1,29) K,VAL(K)                                            
        ENDDO                                                        
      ENDIF
C
      IF(I_AT.EQ.0) WRITE(IUO1,30)                                                          
      IF((IPRINT.GT.0).AND.(I_AT.EQ.0)) THEN
         WRITE(IUO1,31) (IEMET(J),J=1,NEMET) 
      ENDIF  
      ZEM=1.E+20
      DO L=1,NPLAN                                                   
        Z=VAL(L)                                                          
        DO JEMED=1,NEMET                                               
          CALL EMETT(JEMED,IEMET,Z,COORD,NATYP,EMET,NTEM,JNEM,*93)
          IF(I_AT.EQ.0) WRITE(IUO1,34) L,NTEM,EMET(1),EMET(2),EMET(3)
          IF((IPHA.EQ.1).OR.(IPHA.EQ.2)) ZEM=EMET(3) 
          GO TO 33
  93      IF(I_AT.EQ.0) WRITE(IUO1,94) L,NTEM 
  33      CONTINUE
        ENDDO                                          
      ENDDO
C
C..........  Loop on the electrons involved in the  ..........
C..........  spectroscopy :  N_EL = 1 for PHD, XAS  ..........
C..........  LEED or AED and N_EL = 2 for APC       ..........
C 
        DO J_EL=1,N_EL
C
C..........  Writing the information on the spectroscopies  ..........
C..........             in the control file IUO1            ..........
C
        IF(SPECTRO.EQ.'XAS') GOTO 566 
        IF((SPECTRO.EQ.'APC').AND.(J_EL.EQ.1)) THEN
          IF(IPHI.EQ.1) THEN 
            IF(STEREO.EQ.' NO') THEN
              WRITE(IUO1,236)
            ELSE
              WRITE(IUO1,248)
            ENDIF
          ENDIF
          IF(ITHETA.EQ.1) WRITE(IUO1,245)
          IF(I_TEST.EQ.1) WRITE(IUO1,234)
        ENDIF
C
C----------  Photoelectron diffraction case (PHD)  ---------- 
C
        IF((SPECTRO.EQ.'PHD').OR.(SPECTRO.EQ.'APC')) THEN 
          IF(SPECTRO.EQ.'PHD') THEN                 
            IF(IPHI.EQ.1) THEN
              IF(STEREO.EQ.' NO') THEN
                WRITE(IUO1,35)
              ELSE
                WRITE(IUO1,246)
              ENDIF
            ENDIF
            IF(ITHETA.EQ.1) WRITE(IUO1,44) 
            IF(IE.EQ.1) WRITE(IUO1,58)                                         
            IF(INITL.EQ.0) WRITE(IUO1,118)
            IF(I_TEST.EQ.1) WRITE(IUO1,234)
          ENDIF
          IF((SPECTRO.EQ.'APC').AND.(J_EL.EQ.1)) THEN
            WRITE(IUO1,418)
            WRITE(IUO1,18)
          ENDIF
          IF(J_EL.EQ.2) GOTO 222
          IF(IPRINT.GT.0) THEN
            WRITE(IUO1,92)
            WRITE(IUO1,91)
            IF(ISPIN.EQ.0) THEN
              WRITE(IUO1,335)
            ELSE
              WRITE(IUO1,336)
            ENDIF
            WRITE(IUO1,91)
            IF(IPOTC.EQ.0) THEN
              WRITE(IUO1,339)
            ELSE
              WRITE(IUO1,334)
            ENDIF
            WRITE(IUO1,91)
            IF(INITL.NE.0) THEN
              WRITE(IUO1,337)
              WRITE(IUO1,91)
              IF(IPOL.EQ.0) THEN                                                            
                WRITE(IUO1,88)                                                      
              ELSEIF(ABS(IPOL).EQ.1) THEN                                              
                WRITE(IUO1,87)                                                      
              ELSEIF(IPOL.EQ.2) THEN                                              
                WRITE(IUO1,89)                                                      
              ENDIF                                                           
              WRITE(IUO1,91)
              IF(IDICHR.GT.0) THEN
                WRITE(IUO1,338)
              ENDIF 
              WRITE(IUO1,91)
              WRITE(IUO1,92)
              WRITE(IUO1,90)                                                         
              WRITE(IUO1,43) THLUM,PHILUM                            
              IF((SPECTRO.EQ.'PHD').AND.(IMOD.EQ.1)) THEN
                WRITE(IUO1,45)  
              ENDIF                       
            ENDIF
C
            IF(INITL.EQ.2) THEN                                             
              WRITE(IUO1,79) LI,LI-1,LI+1                                      
              IF(I_SO.EQ.1) THEN
                WRITE(IUO1,80) S_O 
              ENDIF 
              DO JE=1,NE                                  
                DO JEM=1,NEMET
                  JTE=IEMET(JEM)
                  IF(ISPIN.EQ.0) THEN
                    WRITE(IUO1,111) JTE,RHOR(JE,JTE,NNL,1,1),
     1                                  RHOR(JE,JTE,NNL,2,1)
                    IF(ITL.EQ.0) THEN
                      WRITE(IUO1,444) JTE,DLT(JE,JTE,NNL,1),
     1                                    DLT(JE,JTE,NNL,2)
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ELSEIF(INITL.EQ.-1) THEN                                        
             WRITE(IUO1,82) LI,LI-1
              IF(I_SO.EQ.1) THEN
                WRITE(IUO1,80) S_O 
              ENDIF 
              DO JE=1,NE                                  
                DO JEM=1,NEMET
                  JTE=IEMET(JEM)
                  IF(ISPIN.EQ.0) THEN
                    WRITE(IUO1,113) JTE,RHOR(JE,JTE,NNL,1,1)
                    IF(ITL.EQ.0) THEN
                      WRITE(IUO1,445) JTE,DLT(JE,JTE,NNL,1)
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ELSEIF(INITL.EQ.1) THEN                                         
              WRITE(IUO1,82) LI,LI+1 
              IF(I_SO.EQ.1) THEN
                WRITE(IUO1,80) S_O 
              ENDIF 
              DO JE=1,NE                                   
                DO JEM=1,NEMET
                  JTE=IEMET(JEM)
                  IF(ISPIN.EQ.0) THEN
                    WRITE(IUO1,113) JTE,RHOR(JE,JTE,NNL,2,1)
                    IF(ITL.EQ.0) THEN
                      WRITE(IUO1,445) JTE,DLT(JE,JTE,NNL,2)
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDIF 
C                                                          
            IF(I_AT.EQ.0) THEN   
              IF(INV(J_EL).EQ.0) THEN                                             
                IF(NDIF.EQ.1) THEN
                  IF(ISPHER.EQ.1) THEN                                      
                    WRITE(IUO1,83)                                                    
                  ELSEIF(ISPHER.EQ.0) THEN                                      
                    WRITE(IUO1,84)                                                    
                  ENDIF                                                         
                ELSE
                  IF(ISPHER.EQ.0) THEN
                    WRITE(IUO1,97) NDIF
                  ELSE
                    WRITE(IUO1,98) NDIF
                  ENDIF  
                ENDIF
              ELSE
                IF(ISPHER.EQ.0) THEN
                  WRITE(IUO1,122)
                ELSE
                  WRITE(IUO1,120) 
                ENDIF  
              ENDIF
            ELSE
              IF(ISPHER.EQ.0) THEN
                WRITE(IUO1,85)                                                    
              ELSE
                WRITE(IUO1,86)                                                    
              ENDIF
            ENDIF
C
          ENDIF
 222      CONTINUE
        ENDIF
C
C----------  LEED case (LED)  ---------- 
C
        IF(SPECTRO.EQ.'LED') THEN 
          IF(IPHI.EQ.1) THEN
            IF(STEREO.EQ.' NO') THEN
              WRITE(IUO1,252)
            ELSE
              WRITE(IUO1,258)
            ENDIF
          ENDIF
          IF(ITHETA.EQ.1) WRITE(IUO1,254) 
          IF(IE.EQ.1) WRITE(IUO1,256)
          IF(IPRINT.GT.0) THEN
            WRITE(IUO1,92)
            WRITE(IUO1,91)
            IF(ISPIN.EQ.0) THEN
              WRITE(IUO1,335)
            ELSE
              WRITE(IUO1,336)
            ENDIF
            WRITE(IUO1,91)
            IF(IPOTC.EQ.0) THEN
              WRITE(IUO1,339)
            ELSE
              WRITE(IUO1,334)
            ENDIF
            WRITE(IUO1,91)
            WRITE(IUO1,92)
            WRITE(IUO1,260)
            WRITE(IUO1,261) THLUM,PHILUM
            IF((SPECTRO.EQ.'LED').AND.(IMOD.EQ.1)) THEN
              WRITE(IUO1,45)
            ENDIF
C
            IF(I_AT.EQ.0) THEN
              IF(INV(J_EL).EQ.0) THEN
                IF(NDIF.EQ.1) THEN
                  IF(ISPHER.EQ.1) THEN
                    WRITE(IUO1,83)
                  ELSEIF(ISPHER.EQ.0) THEN
                    WRITE(IUO1,84)
                  ENDIF
                ELSE
                  IF(ISPHER.EQ.0) THEN
                    WRITE(IUO1,97) NDIF
                  ELSE
                    WRITE(IUO1,98) NDIF
                  ENDIF  
                ENDIF
              ELSE
                IF(ISPHER.EQ.0) THEN
                  WRITE(IUO1,122)
                ELSE
                  WRITE(IUO1,120) 
                ENDIF  
              ENDIF
            ELSE
              IF(ISPHER.EQ.0) THEN
                WRITE(IUO1,85)
              ELSE
                WRITE(IUO1,86)
              ENDIF
            ENDIF
C
          ENDIF
        ENDIF
C
C----------  Auger diffraction case (AED)  ---------- 
C
        IF((SPECTRO.EQ.'AED').OR.(SPECTRO.EQ.'APC')) THEN
          IF(SPECTRO.EQ.'AED') THEN
            IF(IPHI_A.EQ.1) THEN
              IF(STEREO.EQ.' NO') THEN
                WRITE(IUO1,235)
              ELSE
                WRITE(IUO1,247)
              ENDIF
            ENDIF
            IF(ITHETA_A.EQ.1) WRITE(IUO1,244) 
            IF(I_TEST_A.EQ.1) WRITE(IUO1,234)
          ENDIF
          IF((SPECTRO.EQ.'APC').AND.(J_EL.EQ.2)) THEN
            WRITE(IUO1,419)
            WRITE(IUO1,18)
          ENDIF
          IF((SPECTRO.EQ.'AED').OR.(J_EL.EQ.2)) THEN
            IF(IPRINT.GT.0) THEN
              WRITE(IUO1,92)
              WRITE(IUO1,91)
              IF(ISPIN.EQ.0) THEN
                WRITE(IUO1,335)
              ELSE
                WRITE(IUO1,336)
              ENDIF
              WRITE(IUO1,91)
              IF(IPOTC_A.EQ.0) THEN
                WRITE(IUO1,339)
              ELSE
                WRITE(IUO1,334)
              ENDIF
              WRITE(IUO1,91)
              WRITE(IUO1,92)
              WRITE(IUO1,95) AUGER
              CALL AUGER_MULT
              IF(I_MULT.EQ.0) THEN
                WRITE(IUO1,154)
              ELSE
                WRITE(IUO1,155) MULTIPLET
              ENDIF
C
              DO JEM=1,NEMET
                JTE=IEMET(JEM)
                WRITE(IUO1,112) JTE
                DO LE=LE_MIN,LE_MAX
                  WRITE(IUO1,119) LE
                  LA_MIN=L_BOUNDS(LE,1)
                  LA_MAX=L_BOUNDS(LE,2)
                  DO LA=LA_MIN,LA_MAX
                    IF(ISPIN.EQ.0) THEN
                      WRITE(IUO1,115) LA,RHOR_A(LE,JTE,LA,1,1),
     1                                   RHOR_A(LE,JTE,LA,2,1)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
C                                                          
              IF(I_AT.EQ.0) THEN   
                IF(INV(J_EL).EQ.0) THEN                                             
                  IF(NDIF.EQ.1) THEN
                    IF(ISPHER.EQ.1) THEN                                      
                      WRITE(IUO1,83)                                                    
                    ELSEIF(ISPHER.EQ.0) THEN                                      
                      WRITE(IUO1,84)                                                    
                    ENDIF                                                         
                  ELSE
                    IF(ISPHER.EQ.0) THEN
                      WRITE(IUO1,97) NDIF
                    ELSE
                      WRITE(IUO1,98) NDIF
                    ENDIF  
                  ENDIF
                ELSE
                  IF(ISPHER.EQ.0) THEN
                    WRITE(IUO1,122)
                  ELSE
                    WRITE(IUO1,120) 
                  ENDIF  
                ENDIF
              ELSE
                IF(ISPHER.EQ.0) THEN
                  WRITE(IUO1,85)                                                    
                ELSE
                  WRITE(IUO1,86)                                                    
                ENDIF
              ENDIF
C
            ENDIF
          ENDIF
        ENDIF
C
C..........  Check of the dimensioning of the treatment routine  ..........
C
        CALL STOP_TREAT(NFICHLEC,NPLAN,NEMET,NE,NTHETA,NTHETA_A,
     1                  NPHI,NPHI_A,ISOM,I_EXT,I_EXT_A,SPECTRO)
C
C..........     Call of the subroutine performing either     ..........
C..........  the PhD, LEED, AED, EXAFS or APECS calculation  ..........
C
 566    IF(ISPIN.EQ.0) THEN
          IF(SPECTRO.EQ.'PHD') THEN
            CALL PHDDIF_CE(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOR,
     1                     NATCLU,NFICHLEC,JFICH,NP)
          ELSEIF(SPECTRO.EQ.'LED') THEN
c            CALL LEDDIF_CE(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOR,
c     1                     NATCLU,NFICHLEC,JFICH,NP)
          ELSEIF(SPECTRO.EQ.'AED') THEN
c            CALL AEDDIF_CE(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOR_A,
c     1                     NATCLU,NFICHLEC,JFICH,NP,LE_MIN,LE_MAX)
          ELSEIF(SPECTRO.EQ.'XAS') THEN
c            CALL XASDIF_CE(NPLAN,VAL,ZEM,IPHA,RHOR,NFICHLEC,JFICH,NP)
          ELSEIF(SPECTRO.EQ.'APC') THEN
c            IF(J_EL.EQ.1) THEN
c              CALL PHDDIF_CE(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOR,
c     1                    NATCLU,NFICHLEC,JFICH,NP)
c            ELSEIF(J_EL.EQ.2) THEN
c              CALL AEDDIF_CE(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOR_A,
c     1                    NATCLU,NFICHLEC,JFICH,NP,LE_MIN,LE_MAX)
c            ENDIF
          ENDIF
        ELSEIF(ISPIN.EQ.1) THEN
c          IF(SPECTRO.EQ.'PHD') THEN
c            CALL PHDDIF_SP(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOR,
c     1                     NATCLU,NFICHLEC,JFICH,NP)
c          ELSEIF(SPECTRO.EQ.'AED') THEN
c            CALL AEDDIF_SP
c          ELSEIF(SPECTRO.EQ.'XAS') THEN
c            CALL XASDIF_SP
c          ENDIF
         continue
        ENDIF
C
C..........        End of the MS calculation :        ..........
C..........  direct exit or treatment of the results  ..........
C
C
C..........  End of the loop on the electrons  ..........
C
      ENDDO
C
        IF(SPECTRO.EQ.'PHD') THEN
          IF(IPHI.EQ.1) THEN 
            IF(STEREO.EQ.' NO') THEN
              WRITE(IUO1,52)
            ELSE
              WRITE(IUO1,249)
            ENDIF
          ENDIF
          IF(ITHETA.EQ.1) WRITE(IUO1,49)
          IF(IE.EQ.1) WRITE(IUO1,59)
        ELSEIF(SPECTRO.EQ.'LED') THEN
          IF(IPHI.EQ.1) THEN 
            IF(STEREO.EQ.' NO') THEN
              WRITE(IUO1,253)
            ELSE
              WRITE(IUO1,259)
            ENDIF
          ENDIF 
          IF(ITHETA.EQ.1) WRITE(IUO1,255)
          IF(IE.EQ.1) WRITE(IUO1,257)
        ELSEIF(SPECTRO.EQ.'XAS') THEN
          WRITE(IUO1,51) 
        ELSEIF(SPECTRO.EQ.'AED') THEN
          IF(IPHI_A.EQ.1) THEN 
            IF(STEREO.EQ.' NO') THEN
              WRITE(IUO1,237)
            ELSE
              WRITE(IUO1,250)
            ENDIF
          ENDIF 
          IF(ITHETA_A.EQ.1) WRITE(IUO1,238)
        ELSEIF(SPECTRO.EQ.'APC') THEN
          IF(IPHI.EQ.1) THEN 
            IF(STEREO.EQ.' NO') THEN
              WRITE(IUO1,239)
            ELSE
              WRITE(IUO1,251)
            ENDIF
          ENDIF
          IF(ITHETA.EQ.1) WRITE(IUO1,240)
        ENDIF
C
      CLOSE(ICOM)
      IF((NFICHLEC.GT.1).AND.(ISOM.NE.0)) THEN
        WRITE(IUO1,562)
      ENDIF
      IF(ISOM.EQ.0) CLOSE(IUO2)
      IF((ISOM.EQ.0).AND.(NFICHLEC.NE.1)) CLOSE(IUO1)
C
C..........  End of the loop on the data files  ..........
C
      ENDDO
C
      IF(ISOM.NE.0) THEN
        JFF=1
        IF(ISPIN.EQ.0) THEN
          IF(SPECTRO.NE.'XAS') THEN
            CALL TREAT_PHD(ISOM,NFICHLEC,JFF,NP)
          ELSE
c            CALL TREAT_XAS(ISOM,NFICHLEC,NP)
          ENDIF   
        ELSEIF(ISPIN.EQ.1) THEN
c          IF((SPECTRO.EQ.'PHD').OR.(SPECTRO.EQ.'AED')) THEN
c            CALL TREAT_PHD_SP(ISOM,NFICHLEC,JFF,NP)
c          ELSEIF(SPECTRO.EQ.'XAS') THEN
c            CALL TREAT_XAS_SP(ISOM,NFICHLEC,NP)
c          ENDIF
          continue
        ENDIF
      ENDIF
C
      IF((ISOM.NE.0).OR.(NFICHLEC.EQ.1)) CLOSE(IUO1)
      IF(ISOM.NE.0) CLOSE(IUO2)
CST   STOP
      GO TO 999  
C
   1  WRITE(IUO1,60)
      STOP
   2  WRITE(IUO1,61)
      STOP
  55  WRITE(IUO1,65)
      STOP
  56  WRITE(IUO1,64)
      STOP
  74  WRITE(IUO1,75)
      STOP
  99  WRITE(IUO1,100)
      STOP
 180  WRITE(IUO1,181)
      STOP
 182  WRITE(IUO1,183)
      STOP
 184  WRITE(IUO1,185)
      STOP
 504  WRITE(IUO1,505)
      STOP
 510  WRITE(IUO1,511) IUI4
      STOP
 514  WRITE(IUO1,515)
      STOP
 516  WRITE(IUO1,517)
      STOP
 518  WRITE(IUO1,519)
      WRITE(IUO1,889)
      STOP
 520  WRITE(IUO1,521)
      STOP
 540  WRITE(IUO1,541)
      STOP
 550  WRITE(IUO1,551)
      STOP
 570  WRITE(IUO1,571)
      STOP
 580  WRITE(IUO1,581)
      STOP
 590  WRITE(IUO1,591)
      STOP
 600  WRITE(IUO1,601)
      STOP
 602  WRITE(IUO1,603)
      STOP
 604  WRITE(IUO1,605)
      STOP
 606  WRITE(IUO1,607)
      STOP
 608  WRITE(IUO1,609) 
      STOP
 610  WRITE(IUO1,611)
      STOP
 614  WRITE(IUO1,615) NB_AT
      STOP
 620  WRITE(IUO1,621) LE_MAX
      STOP
 630  WRITE(IUO1,631)
      STOP
 890  WRITE(IUO1,891)
      STOP
 895  WRITE(IUO1,896) JA1,JA2
C
   3  FORMAT(5(5X,I4))
   7  FORMAT(3X,F9.4,1X,F9.4,5X,F12.9,5X,F12.9)
CST   9  FORMAT(3X,F9.4,1X,F9.4,5X,E12.6,5X,E12.6)
   9  FORMAT(3X,F9.4,1X,F9.4,E18.6,5X,E18.6)
  17  FORMAT(12X,'ATOM NUMBER ',I4,10X,'CORRESPONDING TRANSLATIONS ',
     1': (',I3,',',I3,',',I3,')')
  18  FORMAT('            ',/)                                    
  20  FORMAT(/,7X,'ATOM OF TYPE ',I2,' AND OF NUMBER ',I5)             
  21  FORMAT(17X,'COORDINATES IN THE TOTAL CLUSTER : (',F7.3,',',
     1       F7.3,',',F7.3,')')                             
  22  FORMAT(22X,'THIS ATOM HAS BEEN SUPRESSED IN THE REDUCED CLUSTER')        
  23  FORMAT(17X,'COORDINATES IN THE REDUCED CLUSTER :(',F7.3,',',
     1       F7.3,',',F7.3,')',5X,'NEW NUMBER : ',I4)                
  24  FORMAT(///,29X,'CONTENTS OF THE REDUCED CLUSTER :',/)            
  26  FORMAT(28X,I4,' ATOMS OF TYPE ',I2)                              
  29  FORMAT(/,20X,'THE Z POSITION OF PLANE ',I3,' IS : ',F6.3)       
  30  FORMAT(///,23X,'THE ABSORBING ATOMS ARE OF TYPE :',/)             
  31  FORMAT(38X,10(I2,3X),//)                                          
  34  FORMAT(//,2X,'PLANE No ',I3,3X,'THE ABSORBER OF TYPE ',  
     1I2,' IS POSITIONED AT (',F7.3,',',F7.3,',',F7.3,')')               
  35  FORMAT(/////,'########## BEGINNING ',      
     1'OF THE AZIMUTHAL PHOTOELECTRON DIFFRACTION CALCULATION #####',   
     2'#####',/////)                                  
  36  FORMAT(/////,'##########              BEGINNING ',      
     1'OF THE EXAFS CALCULATION              ##########',/////)                                  
  37  FORMAT(/////,'++++++++++++++++++++',      
     1'     NUMBERING OF THE ATOMS GENERATED     +++++++++++++++++++')         
  38  FORMAT(///,30X,'TRANSLATION LEVEL : ',I2,///)                 
  39  FORMAT(///,'++++++++++++++++++++++++++++++++++++++++++++++++',    
     1'++++++++++++++++++++++++++++++++',/////)   
  40  FORMAT(/////,'======================',    
     1'     CONTENTS OF THE REDUCED CLUSTER     ======================',
     2      ///)   
  41  FORMAT(///,'====================================================',
     1'============================',/////)
  43  FORMAT(14X,'TH_LIGHT = ',F6.2,' DEGREES',5X,'PHI_LIGHT = ',
     1       F6.2,' DEGREES') 
  44  FORMAT(/////,'########## BEGINNING ',    
     1'OF THE POLAR PHOTOELECTRON DIFFRACTION CALCULATION #####',     
     2'#####',/////)                                  
  45  FORMAT(14X,' (WHEN THE DETECTOR IS ALONG ',
     1       'THE NORMAL TO THE SURFACE)') 
  49  FORMAT(/////,'########## END OF THE ', 
     1'POLAR PHOTOELECTRON DIFFRACTION CALCULATION ##########')                                                
  50  FORMAT(///,22X,'THE CLUSTER IS COMPOSED OF ',I2,' PLANES :')   
  51  FORMAT(/////,'##########                END OF THE ', 
     1'EXAFS CALCULATION                ##########')                                                
  52  FORMAT(/////,'########## END OF THE ', 
     1'AZIMUTHAL PHOTOELECTRON DIFFRACTION CALCULATION #####',
     2'#####')                                             
  57   FORMAT(///,27X,'CALCULATION OF THE SCATTERING FACTOR DONE')        
  58  FORMAT(/////,'########## BEGINNING ',       
     1'OF THE FINE STRUCTURE OSCILLATIONS CALCULATION #####',   
     2'#####',/////)                                  
  59  FORMAT(/////,'########## END OF THE ',   
     1'FINE STRUCTURE OSCILLATIONS CALCULATION #####',
     2'#####')                                                
  60  FORMAT(///,'<<<<<<<<<<  (NAT,NE,NEMET) > (NATP_M,NE_M,',
     1       'NEMET_M) - CHECK THE DIMENSIONING >>>>>>>>>>')
  61  FORMAT(///,22X,' <<<<<<<<<<  THIS STRUCTURE DOES NOT EXIST ',    
     1'  >>>>>>>>>>')                                                   
  64  FORMAT(///,4X,' <<<<<<<<<<  NIV IS TOO SMALL, THE REDUCED ',
     1'CLUSTER HAS NOT CONVERGED YET  >>>>>>>>>>')                          
  65  FORMAT(///,4X,' <<<<<<<<<<  ONLY ONE OF THE VALUES IPHI,ITHETA ', 
     1'ET IE CAN BE EQUAL TO 1  >>>>>>>>>>')                          
  75  FORMAT(///,8X,' <<<<<<<<<<  CHANGE THE DIMENSIONING OF PCREL ',     
     1'IN MAIN ET READ_DATA  >>>>>>>>>>')                                 
  79  FORMAT(//,18X,'INITIAL STATE L = ',I1,5X,'FINAL STATES L = ',      
     1I1,',',I1,/)
  80  FORMAT(15X,'(SPIN-ORBIT COMPONENT OF THE INITIAL CORE STATE : ',
     1          A3,')',//)
  81  FORMAT(18X,'(BOTH SPIN-ORBIT COMPONENTS TAKEN INTO ACCOUNT)')                                                   
  82  FORMAT(//,21X,'INITIAL STATE L = ',I1,5X,'FINAL STATE L = ',I1)     
  83  FORMAT(//,32X,'(SPHERICAL WAVES)')                               
  84  FORMAT(//,34X,'(PLANE WAVES)')                                   
  85  FORMAT(//,26X,'(PLANE WAVES - ATOMIC CASE)')                                  
  86  FORMAT(//,24X,'(SPHERICAL WAVES - ATOMIC CASE)')                                  
  87  FORMAT(24X,'+     LINEARLY POLARIZED LIGHT      +')                                                           
  88  FORMAT(24X,'+        NON POLARIZED LIGHT        +')                                                
  89  FORMAT(24X,'+    CIRCULARLY POLARIZED LIGHT     +') 
  90  FORMAT(////,31X,'POSITION OF THE LIGHT :',/)  
  91  FORMAT(24X,'+',35X,'+') 
  92  FORMAT(24X,'+++++++++++++++++++++++++++++++++++++')                                                       
  94  FORMAT(//,2X,'PLANE No ',I3,3X,'NO ABSORBER OF TYPE ',I2,  
     1' IS PRESENT IN THIS PLANE')
  95  FORMAT(////,31X,'AUGER LINE :',A6,//)
  97  FORMAT(///,19X,'(PLANE WAVES MULTIPLE SCATTERING - ORDER ',I1,
     1       ')')
  98  FORMAT(///,17X,'(SPHERICAL WAVES MULTIPLE SCATTERING - ORDER ',
     1       I1,')')
 100  FORMAT(///,8X,'<<<<<<<<<<  WRONG NAME FOR THE INITIAL STATE',
     1      '  >>>>>>>>>>')
 101  FORMAT(24X,I3,24X,I3)
 102  FORMAT(A1)
 103  FORMAT(31X,F7.2)
 104  FORMAT(29X,F8.5,4X,F8.5,7X,F8.5,4X,F8.5) 
 105  FORMAT(1X,E12.5,1X,E12.5,2X,E12.5,1X,E12.5,4X,E12.5,1X,E12.5,
     1       2X,E12.5,1X,E12.5,2X,E12.5,1X,E12.5,4X,A9)
 106  FORMAT(12X,I3,12X,I3,12X,I3)
 107  FORMAT(5X,I2,5X,I2,5X,I2)
 108  FORMAT(19X,I2,8X,F8.5,1X,F8.5,4X,F8.5,1X,F8.5) 
 109  FORMAT(5X,I2,12X,I2,11X,I2)
 110  FORMAT(16X,'RADIAL MATRIX ELEMENTS FOR THE ABSORBER OF TYPE ',I2,
     1       ' :',/,22X,'(THE SPIN DOUBLET IS GIVEN AS : OUT/IN)',//)
 111  FORMAT(6X,'RADIAL MATRIX ELEMENTS FOR THE ABSORBER OF TYPE ',
     1       I2,' : (',F8.5,',',F8.5,')',/,59X,'(',F8.5,',',F8.5,')')
 112  FORMAT(6X,'RADIAL MATRIX ELEMENTS FOR THE ABSORBER OF TYPE ',
     1       I2,' : ',/,8X,'(LE : ALLOWED VALUES FOR ESCAPING AUGER',
     2      ' ELECTRON)',/,
     2       8X,'(L  : INTERNAL VALUE THAT WILL BE SUMMED ON)',//)
 113  FORMAT(6X,'RADIAL MATRIX ELEMENT FOR THE ABSORBER OF ',
     *       'TYPE  ',I2,' : (',F8.5,',',F8.5,')')
 114  FORMAT(/)
 115  FORMAT(15X,'L = ',I2,5X,'(',F8.5,',',F8.5,')',5X,
     1       '(',F8.5,',',F8.5,')')
 117  FORMAT(12X,I2,5X,I2)
 118  FORMAT(/,37X,'AUGER ELECTRON DIFFRACTION',/)
 119  FORMAT(10X,'LE = ',I2,11X,'DIRECT INTEGRAL',8X,
     1       'EXCHANGE INTEGRAL')
 120  FORMAT(///,15X,'(SPHERICAL WAVES MULTIPLE SCATTERING - MATRIX ',
     1       'INVERSION)')
 122  FORMAT(///,17X,'(PLANE WAVES MULTIPLE SCATTERING - MATRIX ',
     1       'INVERSION)')
 125  FORMAT(11X,A2,5X,I2,3F10.4,12X,I4)
 154  FORMAT(///,20X,'CALCULATION MADE FOR THE FULL AUGER LINE',
     1       ' ',/,' ',/,' ')
 155  FORMAT(///,20X,'CALCULATION MADE FOR THE ',A3,' MULTIPLET ',
     1      'LINE',' ',/,' ',/,' ')
 181  FORMAT(///,'<<<<<<<<<<  NAT OR NE DIFFERENT BETWEEN THE INPUT ',
     1          'AND PHASE SHIFTS FILES  >>>>>>>>>>')
 183  FORMAT(///,'<<<<<<<<<<  NAT OR NE DIFFERENT BETWEEN THE INPUT ',
     1          'AND RADIAL MATRIX ELEMENTS FILES  >>>>>>>>>>')
 185  FORMAT(///,'<<<<<<<<<<  LMAX > NL_M-1 IN THE PHASE SHIFTS ',
     1          'FILE  >>>>>>>>>>')
 234  FORMAT('  ----->  TEST CALCULATION : NO EXCITATION ',
     1           'MATRIX ELEMENTS TAKEN INTO ACCOUNT  <-----',///)
 235  FORMAT(/////,'########## BEGINNING ',      
     1'OF THE AZIMUTHAL AUGER DIFFRACTION CALCULATION #####',   
     2'#####',/////)                                  
 236  FORMAT(/////,'########## BEGINNING ',      
     1'OF THE AZIMUTHAL APECS DIFFRACTION CALCULATION #####',   
     2'#####',/////)                                  
 237  FORMAT(/////,'########## END ',      
     1'OF THE AZIMUTHAL AUGER DIFFRACTION CALCULATION #####',   
     2'#####',/////)                                  
 238  FORMAT(/////,6X,'########## END ',    
     1'OF THE POLAR AUGER DIFFRACTION CALCULATION #####',     
     2'#####',/////)                                  
 239  FORMAT(/////,'########## END ',      
     1'OF THE AZIMUTHAL APECS DIFFRACTION CALCULATION #####',   
     2'#####',/////)                                  
 240  FORMAT(/////,6X,'########## END ',    
     1'OF THE POLAR APECS DIFFRACTION CALCULATION #####',     
     2'#####',/////)                                  
 244  FORMAT(/////,6X,'########## BEGINNING ',    
     1'OF THE POLAR AUGER DIFFRACTION CALCULATION #####',     
     2'#####',/////)                                  
 245  FORMAT(/////,6X,'########## BEGINNING ',    
     1'OF THE POLAR APECS DIFFRACTION CALCULATION #####',     
     2'#####',/////)                                  
 246  FORMAT(/////,'########## BEGINNING ',      
     1'OF THE FULL ANGLE PHOTOELECTRON DIFFRACTION CALCULATION ',
     2'##########',/////)                                  
 247  FORMAT(/////,'########## BEGINNING ',      
     1'OF THE FULL ANGLE AUGER DIFFRACTION CALCULATION ',   
     2'##########',/////)                                  
 248  FORMAT(/////,'########## BEGINNING ',      
     1'OF THE FULL ANGLE APECS DIFFRACTION CALCULATION ',   
     2'##########',/////)                                  
 249  FORMAT(/////,'########## END OF THE ', 
     1'FULL ANGLE PHOTOELECTRON DIFFRACTION CALCULATION #####',
     2'#####')                                             
 250  FORMAT(/////,'########## END ',      
     1'OF THE FULL ANGLE AUGER DIFFRACTION CALCULATION #####',   
     2'#####',/////)                                  
 251  FORMAT(/////,'########## END ',      
     1'OF THE FULL ANGLE APECS DIFFRACTION CALCULATION #####',   
     2'#####',/////)                                  
 252  FORMAT(/////,'########## BEGINNING ',      
     1'OF THE AZIMUTHAL LEED CALCULATION #####',   
     2'#####',/////)                                  
 253  FORMAT(/////,'########## END ',      
     1'OF THE AZIMUTHAL LEED CALCULATION #####',   
     2'#####',/////)                                  
 254  FORMAT(/////,6X,'########## BEGINNING ',    
     1'OF THE POLAR LEED CALCULATION #####',     
     2'#####',/////)                                  
 255  FORMAT(/////,6X,'########## END ',    
     1'OF THE POLAR LEED CALCULATION #####',     
     2'#####',/////)                                  
 256  FORMAT(/////,5X,'########## BEGINNING ',    
     1'OF THE ENERGY LEED CALCULATION #####',     
     2'#####',/////)                                  
 257  FORMAT(/////,5X,'########## END ',    
     1'OF THE ENERGY LEED CALCULATION #####',     
     2'#####',/////)                                  
 258  FORMAT(/////,'########## BEGINNING ',      
     1'OF THE FULL ANGLE LEED CALCULATION ',   
     2'##########',/////)                                  
 259  FORMAT(/////,'########## END OF THE ', 
     1'FULL ANGLE LEED CALCULATION #####',
     2'#####')                                             
 260  FORMAT(////,31X,'POSITION OF THE INITIAL BEAM :',/)  
 261  FORMAT(14X,'TH_BEAM = ',F6.2,' DEGREES',5X,'PHI_BEAM = ',
     1       F6.2,' DEGREES') 
 334  FORMAT(24X,'+   COMPLEX POTENTIAL CALCULATION   +')
 335  FORMAT(24X,'+              STANDARD             +')
 336  FORMAT(24X,'+           SPIN-POLARIZED          +')
 337  FORMAT(24X,'+                WITH               +')
 338  FORMAT(24X,'+         IN DICHROIC MODE          +')
 339  FORMAT(24X,'+    REAL POTENTIAL CALCULATION     +')
 418  FORMAT(///,9X,'------------------------  FIRST ELECTRON : ',
     1       '------------------------')
 419  FORMAT(///,9X,'------------------------ SECOND ELECTRON : ',
     1       '------------------------')
 420  FORMAT(///,9X,'----------------------------------------------',
     1       '----------------------')
 444  FORMAT(12X,'PHASE SHIFTS FOR THE ABSORBER OF TYPE  ',I2,' : ',
     1       '(',F8.5,',',F8.5,')',/,56X,'(',F8.5,',',F8.5,')')
 445  FORMAT(12X,'PHASE SHIFT FOR THE ABSORBER OF TYPE  ',I2,' : (',
     1       F8.5,',',F8.5,')')
 505  FORMAT(///,'<<<<<<<<<<  LI IS LARGER THAN LI_M - ',
     1       'CHECK THE DIMENSIONING  >>>>>>>>>>')
 511  FORMAT(///,'<<<<<<<<<<  NATCLU_M IN THE .inc FILE IS NOT ',
     1       'CONSISTENT WITH THE NUMBER OF ATOMS READ FROM UNIT ',I2,
     2       '  >>>>>>>>>>')
 515  FORMAT(///,'<<<<<<<<<<  INCOMPATIBILITY BETWEEN THE VALUES OF ',
     1       'NAT IN THE DATA AND CLUSTER FILES  >>>>>>>>>>')
 517  FORMAT(///,'<<<<<<<<<<  THERE ARE MISSING VALUES FOR THFWD AND ',
     1       'IBWD  >>>>>>>>>>')
 519  FORMAT(///,'<<<<<<<<<<  NATCLU_M IN THE .inc FILE IS NOT',
     1       ' CONSISTENT WITH THE NUMBER OF ATOMS GENERATED BY THE ',
     2       'CODE  >>>>>>>>>>')
 521  FORMAT(///,'<<<<<<<<<<  SPIN-ORBIT COMPONENT NOT CONSISTENT WITH',
     1       ' THE VALUE OF LI  >>>>>>>>>>')
 530  FORMAT(3X,F9.4,3X,F9.4,3X,F9.4)
 535  FORMAT(29X,F8.5,1X,F8.5) 
 541  FORMAT(///,'<<<<<<<<<<  THE NUMBER OF LINES THFWD DOES NOT ',
     1       'CORRESPOND TO NAT  >>>>>>>>>>') 
 543  FORMAT(5X,F12.9,5X,F12.9)
 549  FORMAT(//,14X,' No ',10X,'COORDINATES',9X,'TYPE',2X,
     2       'SNo',2X,'SYM',/)
 551  FORMAT(///,'<<<<<<<<<<  THE NUMBER OF LINES UJ2 DOES NOT ',
     1       'CORRESPOND TO NAT  >>>>>>>>>>') 
 555  FORMAT(4(7X,I2))
 556  FORMAT(28X,4(I2,5X))
 557  FORMAT(13X,I4,3X,'(',F7.3,',',F7.3,',',F7.3,')',2X,I4,2X,I4,
     1        3X,A2)
 558  FORMAT(/////,18X,'CONTENTS OF THE CLUSTER READ FROM UNIT ',
     1       I2,' : ',/,20X,'READ IN ',A30,//,15X,'No',13X,'(X,Y,Z)',
     2       10X,'CLASS',1X,'ATOM',/)
 559  FORMAT(/////,25X,'CONTENTS OF THE CLUSTER GENERATED : ',//,
     1       14X,' No ',10X,'COORDINATES',9X,'TYPE',2X,'SNo',2X,'SYM',/)
 560  FORMAT(////,12X,'MAXIMAL VALUES OF L FOR THE ',I3,
     1       ' PROTOTYPICAL ATOMS : ',//)
 561  FORMAT(////,18X,'MAXIMAL VALUE OF L FOR THE ',
     1       'PROTOTYPICAL ATOM : ',//)
 562  FORMAT(///,'oooooooooooooooo',12X,'END OF THE INPUT DATA FILE',
     1       13X,'oooooooooooooooo',///)
 563  FORMAT(//,20X,'ENERGY POINT No ',I3,' :',/)
 571  FORMAT(///,'<<<<<<<<<<  THE NUMBER OF LINES ATBAS DOES NOT ',
     1       'CORRESPOND TO NAT  >>>>>>>>>>') 
 581  FORMAT(///,'<<<<<<<<<<  LI OR IMOD NOT CONSISTENT BETWEEN ',
     1       'PHD AND AED FOR COINCIDENCE CALCULATION  >>>>>>>>>>') 
 591  FORMAT(///,'<<<<<<<<<<  THE EXTERNAL DIRECTIONS FILE IS ',
     1       'NOT CONSISTENT WITH THE INPUT DATA FILE  >>>>>>>>>>') 
 601  FORMAT(///,'<<<<<<<<<<  NO_ST_M IS TOO SMALL IN THE .inc FILE ',
     1       '>>>>>>>>>>',//)
 603  FORMAT(///,'<<<<<<<<<<  NSPIN_M OR NSPIN2_M IS TOO SMALL IN THE ',
     1       '.inc FILE  >>>>>>>>>>',//)
 605  FORMAT(///,'<<<<<<<<<<  NT_M IS TOO SMALL IN THE .inc FILE ',
     1       '>>>>>>>>>>',//)
 607  FORMAT(///,'<<<<<<<<<<  THE INITIAL STATE LI IN THE INPUT DATA  ',
     1       'FILE IS DIFFERENT FROM THAT IN THE RADIAL MATRIX ',
     2       'ELEMENTS FILE  >>>>>>>>>>',//)
 609  FORMAT(///,'<<<<<<<<<<  THE TWO TL FILE ARE NOT COMPATIBLE  ',
     1           '>>>>>>>>>>',//)
 611  FORMAT(///,3X,'<<<<<<<<<<  THE RADIAL FILE FOR THE AUGER ',
     1           'ELECTRON IS NOT COMPATIBLE   >>>>>>>>>>',/,
     2           3X,'<<<<<<<<<<  ',17X,'WITH THE INPUT DATA FILE  ',
     3           16X,'>>>>>>>>>>',//)
 613  FORMAT(///,'<<<<<<<<<<  NATP_M SHOULD BE AT LEAST ',I3,'  IN ',
     1       'THE DIMENSIONNING FILE  >>>>>>>>>>',//)
 615  FORMAT(///,'<<<<<<<<<<  NAT_EQ_M SHOULD BE AT LEAST ',I3,'  IN ',
     1       'THE DIMENSIONNING FILE  >>>>>>>>>>',//)
 621  FORMAT(///,'<<<<<<<<<<  LI_M SHOULD BE AT LEAST ',I3,'  IN ',
     1       'THE DIMENSIONNING FILE  >>>>>>>>>>',//)
 631  FORMAT(///,'<<<<<<<<<<      EXCURSIONS OF ANGLES SHOULD ',
     1       ' BE IDENTICAL      >>>>>>>>>>',/,'<<<<<<<<<<     ',
     2       'FOR BOTH ELECTRONS IN CLUSTER ROTATION MODE',
     3        '     >>>>>>>>>>',//)
 776  FORMAT(I2)
 777  FORMAT(A24)
 778  FORMAT(30X,I1)
 779  FORMAT(11X,A2,5X,I2,3F10.4,I5)
 782  FORMAT(/////,22X,'THE CLUSTER GENERATED CONSISTS OF : ',I4,
     1       ' ATOMS')
 889  FORMAT(/////,'<<<<<<<<<<  DECREASE NIV OR INCREASE',
     1       ' NATCLU_M  >>>>>>>>>>')
 891  FORMAT(/////,'<<<<<<<<<<  WRONG NAME FOR THE COORDINATES ''',
     1       'UNITS  >>>>>>>>>>')
 896  FORMAT(///,10X,'<<<<<<<<<<  ERROR IN THE COORDINATES OF THE',
     1       '  ATOMS  >>>>>>>>>>',/,10X,'<<<<<<<<<<     ATOMS ',I4,
     2       ' AND ',I4,' ARE IDENTICAL    >>>>>>>>>>')
C
 999  END                   
