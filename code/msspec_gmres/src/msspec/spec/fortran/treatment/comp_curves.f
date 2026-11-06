C
C
C           ************************************************************
C           * ******************************************************** *
C           * *                                                      * *
C           * *             R-FACTOR, SIMILARITY INDEX,              * *
C           * *              DISTANCE, GOODNESS OF FIT               * * 
C           * *                  KERNEL DISTANCE                     * *
C           * *                 AND SHAPE ANALYSIS                   * *
C           * *                 COMPARISON BETWEEN                   * *
C           * *                AN EXPERIMENTAL FILE                  * *
C           * *             AND RESULTS OF CALCULATIONS              * *
C           * *                                                      * *
C           * ******************************************************** *
C           ************************************************************
C
C  WARNING : the experimental file must contain dots and not commas
C
C
C  Range of the different files :
C
C          Original experimental file       :  1 ---> N_EXP
C          Symmetrized experimental file    :  1 ---> N_SYM 
C          Calculation files                :  1 ---> N_CAL
C          Intersection of files            :  1 ---> N_CMP
C
C
C  Type of comparison:
C
C          I_DIST = 0 ---> R-factor analysis
C          I_DIST = 1 ---> similarity index analysis
C          I_DIST = 2 ---> mathematical distance analysis
C          I_DIST = 3 ---> goodness of fit analysis
C          I_DIST = 4 ---> kernel distance analysis
C
C
C  Name of result files :
C
C          _rfxx.dat  --->  R-factor              (xx from 01 to 12)
C          _sixx.dat  --->  similarity index      (xx from 01 to 12)
C          _dixx.dat  --->  distance              (xx from 01 to 24)
C          _gfxx.dat  --->  goodness of fit       (xx from 01 to 12)
C          _kdxx.dat  --->  kernel distance       (xx from 01 to 12)
C
C
C  Shape analysis : shape descriptors compared instead of points of curve
C
C          I_SA = 0   --->  no shape descriptor used
C                                  I_CUR = 0     raw coordinates used
C                                  I_CUR = 1     modulation function used
C                                  I_CUR = 2     cumulative distribution function used
C                                  I_CUR = 3     curvature function used
C                                  I_CUR = 4     Cui-Femiani-Hu-Wonka-Razdan signature function used
C          I_SA = 1   --->  moments used:
C                                  BASIS = GEOM  geometric moments
C                                  BASIS = LEGE  continuous Legendre 
C                                  BASIS = CHEB  discrete Chebyshev
C                                  BASIS = KRAW  discrete Krawtchouk
C                                  BASIS = HAHN  discrete Hahn
C                                  BASIS = MEIX  discrere Meixner
C                                  BASIS = CHAR  discrete Charlier
C                                  BASIS = SHMA  discrete Shmaliy
C          I_SA = 2   --->  chords used on curve
C                                  I_CHORD = 1   chord_length from point I
C                                  I_CHORD = 2   distance from point I to chord (I-K)-(I+K)
C                                  I_CHORD = 3   chord length along direction THETA
C          I_SA = 3   --->  chain code used on curve 
C                                  N_CONNECT = 3 3-connectivity
C                                  N_CONNECT = 5 5-connectivity
C                                  N_CONNECT = 9 9-connectivity
C          I_SA = 4   --->  curves transformed into closed contour, then:
C                                  SH_AN = CDIS  centroid distance
C                                  SH_AN = TANG  tangent angle
C                                  SH_AN = CURV  curvature function
C                                  SH_AN = TRAR  triangle area
C                                  SH_AN = BEAS  beam angle statistics
C                                  SH_AN = 8CCH  8-connectivity chain code
C                                  SH_AN = CLEN  chord length
C                                  SH_AN = CANG  chord angle
C                                  SH_AN = ACDI  arc chord distance
C                                  SH_AN = FOUR  Fourier descriptors
C
C
C  Note: I_DIST and I_SA can be modified externally through the launching 
C                         script proc_nrfactor
C
C
C  Name of output files :
C
C          _int       --->  interpolated experimental file
C          _sym       --->  symmetrized experimental file       (I_SYM > 0 only)
C          _ren       --->  renormalized/rescaled file
C          _chi       --->  experimental modulation function    (I_CUR = 1 only)
C          _cdf       --->  cumulative distribution function    (I_CUR = 2 or I_DIST = 3 only)
C          _cur       --->  curvature function                  (I_CUR = 3 only)
C          _cfh       --->  Cui-Femiani-Hu-Wonka-Razdan function(I_CUR = 4 only)
C          _wgt       --->  weight function
C          _rec       --->  moment-reconstructed file           (I_SA = 1 only)
C          _mom       --->  moment file                         (I_SA = 1 only)
C          _cco       --->  chain code file                     (I_SA = 2 only)
C          _cho       --->  chord file                          (I_SA = 3 only)
C          _ctr       --->  contour-transformed file            (I_SA = 4 only)
C          _sha       --->  contour shape analysis file         (I_SA = 4 only)
C          _ctr_rec   --->  Fourier reconstructed contour file  (I_SA = 4 only)
C
C
C   Authors :  D. Sébilleau and K. Hatada
C
C  This code is part of the MsSpec package (https://ipr.univ-rennes1.fr/msspec)
C
C  First version: August 2011             Last modified : 19 Jan 2015
C
CST      PROGRAM COMP_CURVE
      SUBROUTINE COMP_CURVES()
C
      PARAMETER (N_SIZE=1000,N_FILES=1000,NMAX=9999)
C
      INTEGER CC_EXP(N_SIZE),CC_CAL(N_SIZE),VALUE,NVALUE
C
      REAL*4 X_CAL(N_SIZE),I_CAL(N_SIZE)
      REAL*4 X_EXP(N_SIZE),I_EXP(N_SIZE)
      REAL*4 I_EXP_1(N_SIZE),I_CAL_1(N_SIZE)
      REAL*4 I_EXP_2(N_SIZE),I_CAL_2(N_SIZE)
      REAL*4 I_EXP_3(N_SIZE),I_CAL_3(N_SIZE)
      REAL*4 I_EXP_4(N_SIZE),I_CAL_4(N_SIZE)
      REAL*4 I_EXP_5(N_SIZE),I_CAL_5(N_SIZE)
      REAL*4 X2_EXP(N_SIZE),I2_EXP(N_SIZE),X2_CAL(N_SIZE),I2_CAL(N_SIZE)
      REAL*4 LARGE,CNORM,EXP_MIN,EXP_MAX
      REAL*4 TEXT1(15),TEXT1B(15),TEXT2(15)
      REAL*4 STEP_EXP,STEP_CAL
      REAL*4 PAR1(N_FILES),PAR2(N_FILES),SYM
      REAL*4 RF1,RF2,RF3,RF4,RF5
      REAL*4 SI1,SI2,SI3,SI4,SI5,SI6,SI7,SI8,SI9,SI10,SI11,SI12
      REAL*4 DI1,DI2,DI3,DI4,DI5,DI6,DI7,DI8,DI9,DI10,DI11,DI12
      REAL*4 DI13,DI14,DI15,DI16,DI17,DI18,DI19,DI20,DI21,DI22,DI23,DI24
      REAL*4 GF1,GF2,GF3,GF4,GF5,GF6,GF7,GF8,GF9,GF10,GF11,GF12
      REAL*4 KD1,KD2,KD3,KD4,KD5,KD6,KD7,KD8,KD9,KD10,KD11,KD12
      REAL*4 M_EXP(0:NMAX),M_CAL(0:NMAX)
      REAL*4 CH_EXP(N_SIZE),CH_CAL(N_SIZE)
      REAL*4 X(N_SIZE),EXPE(N_SIZE),CALC(N_SIZE),W(N_SIZE)
      REAL*4 SIGMA,SHIFT,MAXW 
      REAL*4 CALCULATION(N_SIZE,N_FILES)
C
      REAL*8 MU,NU
C
      CHARACTER*1  STR,CHR,FL
      CHARACTER*3  METHOD,AD
      CHARACTER*4  BASIS,SH_AN,FLAG
      CHARACTER*5  EXT1(12),EXT2(12),EXT3(24),EXT4(12),EXT5(12)
      CHARACTER*40 INFILE,OUTFILE(N_FILES),DUMMY
      CHARACTER*48 RFFILE(24)
      CHARACTER*48 FILE
      CHARACTER*48 CHFILE1,CHFILE2,CHFILE3,CHFILE4,CHFILE5,CHFILE6
      CHARACTER*48 CHFILE7,CHFILE9,CHFILE10,CHFILE11
      CHARACTER*52 CHFILE8
C
      COMMON /PAR_WEI/   I_WEIGHT,I_SHIFT,ALPHA,BETA,SIGMA,MAXW 
C
      DATA ICOM,IUO1        /2,6/
      DATA SMALL,LARGE      /0.0001,1.E+30/
C
      DATA EXT1 /'_rf01','_rf02','_rf03','_rf04','_rf05','_rf06',
     1           '_rf07','_rf08','_rf09','_rf10','_rf11','_rf12'/
      DATA EXT2 /'_si01','_si02','_si03','_si04','_si05','_si06',
     1           '_si07','_si08','_si09','_si10','_si11','_si12'/
      DATA EXT3 /'_di01','_di02','_di03','_di04','_di05','_di06',
     1           '_di07','_di08','_di09','_di10','_di11','_di12',
     2           '_di13','_di14','_di15','_di16','_di17','_di18',
     3           '_di19','_di20','_di21','_di22','_di23','_di24'/
      DATA EXT4 /'_gf01','_gf02','_gf03','_gf04','_gf05','_gf06',
     1           '_gf07','_gf08','_gf09','_gf10','_gf11','_gf12'/
      DATA EXT5 /'_kd01','_kd02','_kd03','_kd04','_kd05','_kd06',
     1           '_kd07','_kd08','_kd09','_kd10','_kd11','_kd12'/
C
C
C..........   Default value for derivatives  ..........
C
C  N_DERIV : number of points used for the calculation 
C               of derivatives (2 <= N_DERIV <= 6)
C
      N_DERIV=3
C
C..........   Initialization of arrays   ..........
C
      DO J=1,N_SIZE
C
        X_EXP(J)=0.
        I_EXP(J)=0.
        X_CAL(J)=0.
        I_CAL(J)=0.
C
        I_EXP_1(J)=0.
        I_EXP_2(J)=0.
        I_EXP_3(J)=0.
        I_EXP_4(J)=0.
        I_EXP_5(J)=0.
        I_CAL_1(J)=0.
        I_CAL_2(J)=0.
        I_CAL_3(J)=0.
        I_CAL_4(J)=0.
        I_CAL_5(J)=0.
C
        X2_EXP(J)=0.
        I2_EXP(J)=0.
        X2_CAL(J)=0.
        I2_CAL(J)=0.
C
        CH_EXP(J)=0.
        CH_CAL(J)=0.
C
        X(J)=0.
        EXPE(J)=0.
        CALC(J)=0.
C
        CC_EXP(J)=0
        CC_CAL(J)=0
C
      ENDDO
C
      DO J=1,N_FILES
C
        PAR1(J)=0.
        PAR2(J)=0.
C
      ENDDO
C
      DO J=0,NMAX
C
        M_EXP(J)=0.
        M_CAL(J)=0.
C
      ENDDO
C
C..........   Opening the input data file
C
      OPEN(UNIT=ICOM, FILE='comp_curves.dat', STATUS='OLD')
C
C  Checking for external flag in input data file that will 
C    trigger reading of external data in the script file
C
      READ(ICOM,51) FL
      BACKSPACE ICOM
C
      IF(FL.EQ.'+') THEN
        I_SCRI=1
      ELSE
        I_SCRI=0
      ENDIF
C
C  Reading external parameters from the launching script file
C
      IF(I_SCRI.EQ.1) THEN
        READ(*,*,ERR=6) FLAG
        READ(*,*,ERR=6) NVALUE
        IF(FLAG.EQ.'DIST') THEN
          I_EXT=1
          CHR=CHAR(NVALUE+48)
          AD='_d'//CHR
        ELSEIF(FLAG.EQ.'SHAP') THEN
          I_EXT=1
          CHR=CHAR(NVALUE+48)
          AD='_s'//CHR
        ELSE
          I_EXT=0
        ENDIF
C
  6     CONTINUE
      ELSE
        I_EXT=0
      ENDIF
C
C..........   Opening the output check file
C
CST      IF(I_EXT.EQ.0) THEN
CST        OPEN(UNIT=IUO1, FILE='comp_curves.lis', STATUS='UNKNOWN')
CST      ELSE
CST        OPEN(UNIT=IUO1, FILE='comp_curves.lis', ACCESS='APPEND',
CST     1       STATUS='UNKNOWN')
CST      ENDIF
C
      IF(I_SCRI.EQ.1) THEN
        WRITE(IUO1,464)
      ENDIF
C
C..........   Reading of the input data   ..........
C
C
      READ(ICOM,1) DUMMY
      READ(ICOM,2) TEXT1
      READ(ICOM,2) TEXT1B
      READ(ICOM,1) DUMMY
C
C         General parameters
C
      READ(ICOM,1) DUMMY
      READ(ICOM,2) TEXT2
      READ(ICOM,1) DUMMY
C
      READ(ICOM,14) N_PAR,NORM,I_SCALE,I_NORM
      READ(ICOM,15) I_SYM,SYM,I_POSI
      READ(ICOM,16) I_DIST,I_CUR,I_SA,I_PRINT
C
C         Taking into account external changes
C
      IF(I_EXT.EQ.1) THEN
        IF(FLAG.EQ.'DIST') THEN
          I_DIST=NVALUE
          WRITE(IUO1,462) NVALUE
        ELSEIF(FLAG.EQ.'SHAP') THEN
          I_SA=NVALUE
          WRITE(IUO1,463) NVALUE
        ENDIF
      ENDIF
C
C         Weight parameters
C
      DO LINE=1,3
        READ(ICOM,1) DUMMY
      ENDDO
C
      READ(ICOM,17) I_WEIGHT,ALPHA,BETA,SIGMA
      READ(ICOM,18) I_SHIFT,MAXW
C
C         R-factors parameter
C
      DO LINE=1,3
        READ(ICOM,1) DUMMY
      ENDDO
C
      READ(ICOM,19) V_I
C
C         Similarity indices parameters
C
      DO LINE=1,3
        READ(ICOM,1) DUMMY
      ENDDO
C
      READ(ICOM,20) ALPHAS,BETAS,N_BINS
C
C         Distances parameters
C
      DO LINE=1,3
        READ(ICOM,1) DUMMY
      ENDDO
C
      READ(ICOM,21) ALPHAD,I_BETA,L,SIGMAD 
C
C         Goodness of fit parameters
C
      DO LINE=1,3
        READ(ICOM,1) DUMMY
      ENDDO
C
      READ(ICOM,29) N_BING,ALPHAG
C
C         Kernel distances parameters
C
      DO LINE=1,3
        READ(ICOM,1) DUMMY
      ENDDO
C
      READ(ICOM,22) ALPHAK,L,SIGMAK
C
C         Moments parameters
C
      DO LINE=1,3
        READ(ICOM,1) DUMMY
      ENDDO
C
      READ(ICOM,23) N_GRID,N_MOM,BASIS  
      READ(ICOM,24) I_ALG,MU,NU
C
C         Chords parameters
C
      DO LINE=1,3
        READ(ICOM,1) DUMMY
      ENDDO
C
      READ(ICOM,25) I_CHORD,METHOD,VALUE,N_BINC
C
C         Chain codes parameter
C
      DO LINE=1,3
        READ(ICOM,1) DUMMY
      ENDDO
C
      READ(ICOM,26) N_CONNECT,SCALEC
C
C         Contour parameters
C
      DO LINE=1,3
        READ(ICOM,1) DUMMY
      ENDDO
C
      READ(ICOM,27) NBIN,N_LEN,SH_AN,I_FOU
      READ(ICOM,28) INORM
C
      DO LINE=1,5
        READ(ICOM,1) DUMMY
      ENDDO
C
C         Experimental file
C
      READ(ICOM,34) INFILE
C
      DO LINE=1,5
        READ(ICOM,1) DUMMY
      ENDDO
C
C..........   Checking the number of calculations   .........
C
      NFILE=0
      DO JLINE=1,N_FILES
        READ(ICOM,3) STR
        IF(STR.EQ.'+') THEN
           GOTO 5
        ELSE
          NFILE=NFILE+1
        ENDIF
      ENDDO
C
      IF(NFILE.GT.N_FILES) THEN
        WRITE(IUO1,11) NFILE
        STOP
      ENDIF
C
   5  REWIND ICOM
      DO LINE=1,64
        READ(ICOM,1) DUMMY
      ENDDO
C
      DO LINE=1,NFILE
        READ(ICOM,35) OUTFILE(LINE),PAR1(LINE),PAR2(LINE)
      ENDDO
C
      CLOSE(ICOM)
C
C
C..........   Writing the input data   ..........
C..........     into the check file    ..........
C
C
      WRITE(IUO1,100)
      WRITE(IUO1,101)
      WRITE(IUO1,101)
      WRITE(IUO1,102) TEXT1
      WRITE(IUO1,102) TEXT1B
      WRITE(IUO1,101)
      WRITE(IUO1,101)
      WRITE(IUO1,203)
C
      WRITE(IUO1,140)
      WRITE(IUO1,114) N_PAR,NORM,I_SCALE,I_NORM
      WRITE(IUO1,115) I_SYM,SYM,I_POSI
      WRITE(IUO1,116) I_DIST,I_CUR,I_SA,I_PRINT
C
      WRITE(IUO1,141)
      WRITE(IUO1,117) I_WEIGHT,ALPHA,BETA,SIGMA
      WRITE(IUO1,118) I_SHIFT,MAXW
C
      IF(I_DIST.EQ.0) THEN
        WRITE(IUO1,142)
        WRITE(IUO1,119)  V_I
      ENDIF
C
      IF(I_DIST.EQ.1) THEN
        WRITE(IUO1,143)
        WRITE(IUO1,120)  ALPHAS,BETAS,N_BINS
      ENDIF
C
      IF(I_DIST.EQ.2) THEN
        WRITE(IUO1,144)
        WRITE(IUO1,121) ALPHAD,I_BETA,L,SIGMAD
      ENDIF
C
      IF(I_DIST.EQ.3) THEN
        WRITE(IUO1,155)
        WRITE(IUO1,129) N_BING,ALPHAG
      ENDIF
C
      IF(I_DIST.EQ.4) THEN
        WRITE(IUO1,145)
        WRITE(IUO1,122) ALPHAK,L,SIGMAK
      ENDIF
C
      IF(I_SA.EQ.1) THEN
        WRITE(IUO1,146)
        WRITE(IUO1,123) N_GRID,N_MOM,BASIS  
        WRITE(IUO1,124) I_ALG,MU,NU
      ENDIF
C
      IF(I_SA.EQ.2) THEN
        WRITE(IUO1,147)
        WRITE(IUO1,125) I_CHORD,METHOD,VALUE,N_BINC
      ENDIF
C
      IF(I_SA.EQ.3) THEN
        WRITE(IUO1,148)
        WRITE(IUO1,126) N_CONNECT,SCALEC
      ENDIF
C
      IF(I_SA.EQ.4) THEN
        WRITE(IUO1,149)
        WRITE(IUO1,127) NBIN,N_LEN,SH_AN,I_FOU
        WRITE(IUO1,128) INORM
      ENDIF
C
C..........   Writing the type of normalization   ..........
C..........   and the type of analysis performed  ..........
C
      IF(I_DIST.EQ.0) THEN
        WRITE(IUO1,400)
      ELSEIF(I_DIST.EQ.1) THEN 
        WRITE(IUO1,401)
      ELSEIF(I_DIST.EQ.2) THEN 
        WRITE(IUO1,402)
      ELSEIF(I_DIST.EQ.3) THEN 
        IF(I_POSI.EQ.1) THEN
          WRITE(IUO1,403)
        ELSEIF(I_POSI.EQ.0) THEN
           WRITE(IUO1,458)
           I_POSI=1
        ENDIF
      ELSEIF(I_DIST.EQ.4) THEN 
        WRITE(IUO1,404)
      ENDIF
C
      IF(I_SA.EQ.0) THEN
        IF(I_CUR.EQ.0) THEN
          WRITE(IUO1,406)
        ELSEIF(I_CUR.EQ.1) THEN
          WRITE(IUO1,405)
        ELSEIF(I_CUR.EQ.2) THEN
          WRITE(IUO1,421)
          IF(NORM.GT.0) THEN
            WRITE(IUO1,459)
            NORM=0
          ENDIF
          IF(I_SCALE.EQ.0) THEN
            WRITE(IUO1,460)
            I_SCALE=1
          ENDIF
        ELSEIF(I_CUR.EQ.3) THEN
          WRITE(IUO1,422)
        ELSEIF(I_CUR.EQ.4) THEN
          WRITE(IUO1,423)
        ENDIF
      ELSEIF(I_SA.EQ.1) THEN
        IF(BASIS.EQ.'GEOM') WRITE(IUO1,407)
        IF(BASIS.EQ.'LEGE') WRITE(IUO1,408)
        IF(BASIS.EQ.'CHEB') WRITE(IUO1,409)
        IF(BASIS.EQ.'KRAW') WRITE(IUO1,410)
        IF(BASIS.EQ.'HAHN') WRITE(IUO1,411)
        IF(BASIS.EQ.'MEIX') WRITE(IUO1,412)
        IF(BASIS.EQ.'CHAR') WRITE(IUO1,413)
        IF(BASIS.EQ.'SHMA') WRITE(IUO1,414)
      ELSEIF(I_SA.EQ.2) THEN
         IF(I_CHORD.EQ.1) WRITE(IUO1,415)
         IF(I_CHORD.EQ.2) WRITE(IUO1,416)
         IF(I_CHORD.EQ.3) WRITE(IUO1,417)
      ELSEIF(I_SA.EQ.3) THEN
         IF(N_CONNECT.EQ.3) WRITE(IUO1,418) 
         IF(N_CONNECT.EQ.5) WRITE(IUO1,419) 
         IF(N_CONNECT.EQ.9) WRITE(IUO1,420) 
      ELSEIF(I_SA.EQ.4) THEN
        IF(SH_AN.EQ.'CDIS') WRITE(IUO1,437)
        IF(SH_AN.EQ.'TANG') WRITE(IUO1,438)
        IF(SH_AN.EQ.'CURV') WRITE(IUO1,439)
        IF(SH_AN.EQ.'TRAR') WRITE(IUO1,440)
        IF(SH_AN.EQ.'BEAS') WRITE(IUO1,441)
        IF(SH_AN.EQ.'8CCH') WRITE(IUO1,442)
        IF(SH_AN.EQ.'CLEN') WRITE(IUO1,443)
        IF(SH_AN.EQ.'CANG') WRITE(IUO1,444)
        IF(SH_AN.EQ.'ACDI') WRITE(IUO1,445)
        IF(SH_AN.EQ.'FOUR') WRITE(IUO1,446)
        IF(NORM.GT.0) THEN
          WRITE(IUO1,461)
          NORM=0
        ENDIF
      ENDIF
C
      IF((I_SA.GT.0).AND.(I_NORM.GT.0)) THEN
        WRITE(IUO1,457)
        I_NORM=0
      ENDIF
C
      IF(NORM.EQ.0) THEN
        WRITE(IUO1,450)
      ELSEIF(NORM.EQ.1) THEN
        WRITE(IUO1,451)
      ELSEIF(NORM.EQ.2) THEN
        WRITE(IUO1,452)
      ELSEIF(NORM.EQ.3) THEN
        WRITE(IUO1,453)
      ELSEIF(NORM.EQ.4) THEN
        WRITE(IUO1,454)
      ELSEIF(NORM.EQ.5) THEN
        WRITE(IUO1,455)
      ENDIF
C
      IF(I_SCALE.EQ.1) THEN
        WRITE(IUO1,456)
      ENDIF
C
C..........   Writing the names of experimental   ..........
C..........         and calculation files         ..........
C..........        against the parameters         ..........
C
      IF(N_PAR.EQ.1) THEN
        WRITE(IUO1,150) INFILE
        WRITE(IUO1,151) OUTFILE(1),PAR1(1)
        DO JFILE=2,NFILE
          WRITE(IUO1,152) OUTFILE(JFILE),PAR1(JFILE)
        ENDDO 
      ELSEIF(N_PAR.EQ.2) THEN
        WRITE(IUO1,153) INFILE
        WRITE(IUO1,154) OUTFILE(1),PAR1(1),PAR2(1)
        DO JFILE=2,NFILE
          WRITE(IUO1,152) OUTFILE(JFILE),PAR1(JFILE),PAR2(JFILE)
        ENDDO 
      ENDIF
C
C..........   Checking the number of moments   ......... 
C
      IF(I_SA.EQ.1) THEN
        IF(N_MOM.GT.NMAX) THEN
          WRITE(IUO1,308) 
          STOP
        ENDIF
      ENDIF
C
C..........   Checking the consistency of the calculated files   ......... 
C..........     (must have same lower, upper bounds and step)    .........
C
      IF(NFILE.GT.1) CALL CHECK_CALC_FILE(OUTFILE,NFILE,IUO1,STEP_CAL)
C
C..........   Number of R-factors/Similarity/etc to compute   ......... 
C
      IF(I_DIST.EQ.0) THEN
        N_RF=12
      ELSEIF(I_DIST.EQ.1) THEN
        N_RF=12
      ELSEIF(I_DIST.EQ.2) THEN
        N_RF=24
      ELSEIF(I_DIST.EQ.3) THEN
        N_RF=12
      ELSEIF(I_DIST.EQ.4) THEN
        N_RF=12
      ENDIF
C
C  Finding the real size of the experimental file name 
C              and the position of the dot
C
      N_DOT=1
      DO J_CHAR=1,40
        IF(INFILE(J_CHAR:J_CHAR).EQ.'.') GOTO 50
        N_DOT=N_DOT+1
      ENDDO
  50  CONTINUE
C
      N_CHAR=0
      DO J_CHAR=1,40
        IF(INFILE(J_CHAR:J_CHAR).EQ.' ') GOTO 500
        N_CHAR=N_CHAR+1
      ENDDO
 500  CONTINUE
C
C..........   Opening the R-factor/Similarity/   ........
C..........     Distance/Goodness of fit/        .........
C..........       Kernel distance files          .........
C
      DO JR=1,N_RF
        IF(I_DIST.EQ.0) THEN
          IF(I_EXT.EQ.0) THEN
            RFFILE(JR)='rfc'//INFILE(4:N_DOT-1)//EXT1(JR)//
     1                        INFILE(N_DOT:N_CHAR)
          ELSE
            RFFILE(JR)='rfc'//INFILE(4:N_DOT-1)//EXT1(JR)//AD//
     1                        INFILE(N_DOT:N_CHAR)
          ENDIF
        ELSEIF(I_DIST.EQ.1) THEN
          IF(I_EXT.EQ.0) THEN
            RFFILE(JR)='rfc'//INFILE(4:N_DOT-1)//EXT2(JR)//
     1                        INFILE(N_DOT:N_CHAR)
          ELSE
            RFFILE(JR)='rfc'//INFILE(4:N_DOT-1)//EXT2(JR)//AD//
     1                        INFILE(N_DOT:N_CHAR)
          ENDIF
        ELSEIF(I_DIST.EQ.2) THEN
          IF(I_EXT.EQ.0) THEN
            RFFILE(JR)='rfc'//INFILE(4:N_DOT-1)//EXT3(JR)//
     1                        INFILE(N_DOT:N_CHAR)
          ELSE
            RFFILE(JR)='rfc'//INFILE(4:N_DOT-1)//EXT3(JR)//AD//
     1                        INFILE(N_DOT:N_CHAR)
          ENDIF
        ELSEIF(I_DIST.EQ.3) THEN
          IF(I_EXT.EQ.0) THEN
            RFFILE(JR)='rfc'//INFILE(4:N_DOT-1)//EXT4(JR)//
     1                        INFILE(N_DOT:N_CHAR)
          ELSE
            RFFILE(JR)='rfc'//INFILE(4:N_DOT-1)//EXT4(JR)//AD//
     1                        INFILE(N_DOT:N_CHAR)
          ENDIF
        ELSEIF(I_DIST.EQ.4) THEN
          IF(I_EXT.EQ.0) THEN
            RFFILE(JR)='rfc'//INFILE(4:N_DOT-1)//EXT5(JR)//
     1                        INFILE(N_DOT:N_CHAR)
          ELSE
            RFFILE(JR)='rfc'//INFILE(4:N_DOT-1)//EXT5(JR)//AD//
     1                        INFILE(N_DOT:N_CHAR)
          ENDIF
        ENDIF
        NUNIT=JR+10
        OPEN(UNIT=NUNIT, FILE=RFFILE(JR), STATUS='unknown')
      ENDDO
      NU_LAST=NUNIT
C
C  Names of the interpolated/symmetrized experimental files,
C    names of the transformed files (moments, contours, ...),
C    and name of R-factor/Similarity/Distance/Goodness of fit/Kernel files
C
      NUNIT3=NU_LAST+3
      NUNIT4=NU_LAST+4
      NUNIT5=NU_LAST+5
      NUNIT6=NU_LAST+6
      NUNIT7=NU_LAST+7
      NUNIT9=NU_LAST+9
      NUNIT10=NU_LAST+10
C
      I_SW=1
      CALL NAME_FILES(I_SW,I_DIST,I_SYM,I_CUR,I_SA,NU_LAST,SH_AN,INFILE,
     1                CHFILE1,CHFILE2,CHFILE3,CHFILE4,CHFILE5,CHFILE6,
     2                CHFILE7,CHFILE8,CHFILE9,CHFILE10,CHFILE11)
C
C..........  Intersection between experimental and calculation files   .........
C
      CALL FILE_INTERSECTION(NU_LAST,INFILE,OUTFILE(1),IUO1,
     1                       N_EXP,N_CAL,STEP_EXP,STEP_CAL,N_CMP,J_IN,
     2                       J_FI,X_EXP,I_EXP,X_CAL,I_CAL)
C
C  Symmetrization of experimental file whenever necessary
C
      IF(I_SYM.NE.0) THEN
        CALL SYMMETRIZE(X_EXP,I_EXP,N_EXP,I_SYM,SYM,STEP_EXP,IUO1,
     1                  NUNIT3,CHFILE3)
      ENDIF
C
C  Interpolation of the (symmetrized) experimental file on the 
C    on the calculation grid within the intersection bounds.
C
C         The result is stored as (X2_EXP,I2_EXP)
C
      CALL EXPE_INTERPOLATE(X_EXP,I_EXP,X_CAL,I_CAL,N_EXP,N_CMP,J_IN,
     1                      STEP_CAL,IUO1,NUNIT3,X2_EXP,I2_EXP)
C
C   Normalization of the experimental curve whenever required
C
      IF((NORM.GT.0).AND.(I_CUR.NE.1)) THEN
        CALL NORMALIZE_CURVE(X2_EXP,I2_EXP,N_CMP,NORM,IUO1)
      ENDIF
C
C..............  Shape analysis or not on experimental file:  ..............
C..............      selection of WHAT is to be compared      ..............
C
      IF(I_SA.EQ.0) THEN
C
C   Computing the experimental modulation function
C
        IF(I_CUR.EQ.1) THEN
          CALL MODULATION_FUNCTION(X2_EXP,I2_EXP,N_CMP,NUNIT6,IUO1,
     1                             NORM,0,I_PRINT)
C
C   Calculating the experimental cumulative distribution function
C
        ELSEIF(I_CUR.EQ.2) THEN
          CALL DISTRIBUTION_FUNCTION(X2_EXP,I2_EXP,N_CMP,NUNIT5)
C
C   Calculating the experimental curvature function
C
        ELSEIF(I_CUR.EQ.3) THEN
          CALL CURVATURE_FUNCTION(X2_EXP,I2_EXP,N_CMP,NUNIT9)
C
C   Calculating the experimental Cui-Femiani-Hu-Wonka-Razdan signature function
C
        ELSEIF(I_CUR.EQ.4) THEN
          CALL CFHWR_FUNCTION(X2_EXP,I2_EXP,N_CMP,NUNIT10)
C
        ENDIF
C
C   Computing the moments of the experimental file
C
      ELSEIF(I_SA.EQ.1) THEN
C
        CALL MOMENTS(I2_EXP,X2_EXP,N_CMP,N_GRID,IUO1,I_ALG,MU,NU,BASIS,
     1               M_EXP,N_MOM)
        WRITE(IUO1,49) CHFILE7
C
C   Computing the chords of the experimental file
C
      ELSEIF(I_SA.EQ.2) THEN
C
        CALL CHORDS(X2_EXP,I2_EXP,N_CMP,I_CHORD,METHOD,VALUE,N_BINC,
     1              IUO1,CH_EXP)
        WRITE(IUO1,48) CHFILE7
C
C   Computing the chain codes of the experimental file
C
      ELSEIF(I_SA.EQ.3) THEN 
C
        CALL CHAIN_CODE(X2_EXP,I2_EXP,N_CMP,N_CONNECT,SCALEC,IUO1,
     1                  CC_EXP) 
        WRITE(IUO1,47) CHFILE7
C
C   Transforming the experimental file into a closed contour
C     in order to use contour analysis methods
C
      ELSEIF(I_SA.EQ.4) THEN
C
C......  Checking concistency of NBIN and N_LEN input parameters
C......      with respect of the number of grid points N_CMP 
C
        IF(SH_AN.EQ.'TRAR') THEN
          NMX=N_CMP/2-1
          IF(N_LEN.GT.NMX) THEN
            WRITE(IUO1,309) NMX
            STOP
          ENDIF
        ELSEIF(SH_AN.EQ.'BEAS') THEN
          NMX=N_CMP/2-1
          IF(N_LEN.GT.NMX) THEN
            WRITE(IUO1,309) NMX
            STOP
          ENDIF
        ELSEIF(SH_AN.EQ.'CLEN') THEN
          NMX=N_CMP-1
          IF(N_LEN.EQ.0) THEN
            IF(NBIN.GT.NMX) THEN
              WRITE(IUO1,310) NMX
              STOP
            ENDIF
          ELSE
            IF(N_LEN.GT.NMX) THEN
              WRITE(IUO1,309) NMX
              STOP
            ENDIF
          ENDIF
        ELSEIF(SH_AN.EQ.'CANG') THEN
          NMX=N_CMP-1
          IF(N_LEN.EQ.0) THEN
            IF(NBIN.GT.NMX) THEN
              WRITE(IUO1,310) NMX
              STOP
            ENDIF
          ELSE
            IF(N_LEN.GT.NMX) THEN
              WRITE(IUO1,309) NMX
              STOP
            ENDIF
          ENDIF
        ELSEIF(SH_AN.EQ.'ACDI') THEN
          NMX=N_CMP-1
          IF(N_LEN.GT.NMX) THEN
            WRITE(IUO1,309) NMX
            STOP
          ENDIF
        ENDIF
C
        FILE=CHFILE6
C
        CALL CONTOUR(I2_EXP,X2_EXP,N_CMP,IUO1,NBIN,N_LEN,SH_AN,
     1               I_FOU,INORM,FILE,M_EXP)
        WRITE(IUO1,36) CHFILE7
C
        CLOSE(99)
        CLOSE(98)
        CLOSE(97)
C
      ENDIF
C
C..............  End of shape analysis on experimental file  ..............
C
C
C    Opening the calculation files and printing the step
C
      IF(I_PRINT.EQ.1) WRITE(IUO1,312) STEP_CAL
C
C  Calculation of the maximum and minimum of experiment
C
      EXP_MIN=LARGE
      EXP_MAX=0.
C
      DO J=1,N_CMP
        EXP_MIN=MIN(EXP_MIN,I2_EXP(J))
        EXP_MAX=MAX(EXP_MAX,I2_EXP(J))
      ENDDO
C
C  Opening of the calculated files and storing their descriptors
C             in the CALCULATION array
C
      NUNIT3=NU_LAST+3
      NUNIT4=NU_LAST+4
      NUNIT5=NU_LAST+5
      NUNIT6=NU_LAST+6
      NUNIT7=NU_LAST+7
      NUNIT9=NU_LAST+9
      NUNIT10=NU_LAST+10
C
      DO JFILE=1,NFILE
C
        OPEN(UNIT=NUNIT3, FILE=OUTFILE(JFILE), STATUS='unknown')
C
C  Names of the transformed files for shape analysis
C
        I_SW=2
        CALL NAME_FILES(I_SW,I_DIST,I_SYM,I_CUR,I_SA,NU_LAST,SH_AN,
     1                  OUTFILE(JFILE),CHFILE1,CHFILE2,CHFILE3,CHFILE4,
     2                  CHFILE5,CHFILE6,CHFILE7,CHFILE8,CHFILE9,
     3                  CHFILE10,CHFILE11)
C
C  Calculation of the maximum, minimum and mean
C
         DO JLINE=1,N_CAL
           READ(NUNIT3,*) X_CAL(JLINE),I_CAL(JLINE)
         ENDDO
C
        CLOSE(NUNIT3)
C
C  Reshuffling the x values of the calculation file
C    to keep only those of the intersection
C
        JJ=0
        DO J=J_IN,J_FI
          JJ=JJ+1
          X2_CAL(JJ)=X_CAL(J)
          I2_CAL(JJ)=I_CAL(J)
        ENDDO
C
        IF((NORM.GT.0).AND.(I_CUR.NE.1)) THEN
C
C   Normalization of the calculated curve whenever required
C
          CALL NORMALIZE_CURVE(X2_CAL,I2_CAL,N_CMP,NORM,IUO1)
C
        ENDIF
C
        IF(I_SCALE.EQ.1) THEN
C
C  Scaling of the calculation to the min and max of the experiment
C
          CALL RESCALE_TO_EXP(I2_CAL,N_CMP,EXP_MIN,EXP_MAX)
C
        ENDIF
C
C..............  Shape analysis or not on calculated file:  ..............
C..............     selection of WHAT is to be compared     ..............
C
        IF(I_SA.EQ.0) THEN
C
C  Calculation of the modulation function
C      I2_CAL becomes this function
C
          IF(I_CUR.EQ.1) THEN
C
            CALL MODULATION_FUNCTION(X2_CAL,I2_CAL,N_CMP,NUNIT6,IUO1,
     1                               NORM,1,I_PRINT)
C
C   Computing the calculation cumulative distribution function
C      I2_CAL becomes this function
C
          ELSEIF(I_CUR.EQ.2) THEN
            CALL DISTRIBUTION_FUNCTION(X2_CAL,I2_CAL,N_CMP,NUNIT5)
C
C   Computing the calculation curvature function
C      I2_CAL becomes this function
C
          ELSEIF(I_CUR.EQ.3) THEN
            CALL CURVATURE_FUNCTION(X2_CAL,I2_CAL,N_CMP,NUNIT9)
C
C   Computing the calculation Cui-Femiani-Hu-Wonka-Razdan signature function
C      I2_CAL becomes this function
C
          ELSEIF(I_CUR.EQ.4) THEN
            CALL CFHWR_FUNCTION(X2_CAL,I2_CAL,N_CMP,NUNIT10)
C
          ENDIF
C
        ELSEIF(I_SA.EQ.1) THEN
C
C   Computing the moments of the calculated file
C
          CALL MOMENTS(I2_CAL,X2_CAL,N_CMP,N_GRID,IUO1,I_ALG,MU,NU,
     1                 BASIS,M_CAL,N_MOM)
          WRITE(IUO1,49) CHFILE7
C    
        ELSEIF(I_SA.EQ.2) THEN
C
C   Computing the chords of the calculated file
C
          CALL CHORDS(X2_CAL,I2_CAL,N_CMP,I_CHORD,METHOD,VALUE,N_BINC,
     1                IUO1,CH_CAL)
          WRITE(IUO1,48) CHFILE7
C    
        ELSEIF(I_SA.EQ.3) THEN
C
C   Computing the chain codes of the calculated file
C
          CALL CHAIN_CODE(X2_CAL,I2_CAL,N_CMP,N_CONNECT,SCALEC,IUO1,
     1                    CC_CAL) 
          WRITE(IUO1,47) CHFILE7
C    
        ELSEIF(I_SA.EQ.4) THEN
C
C   Transforming the calculated file into a closed contour
C     in order to use contour analysis methods
C
          FILE=CHFILE6
C
          CALL CONTOUR(I2_CAL,X2_CAL,N_CMP,IUO1,NBIN,N_LEN,SH_AN,
     1                 I_FOU,INORM,FILE,M_CAL)
          WRITE(IUO1,36) CHFILE7
C
          CLOSE(99)
          CLOSE(98)
          CLOSE(97)
C
        ENDIF
C
C..............  End of shape analysis on calculated file  ..............
C
C
C    Storage of the data to be compared into X,EXPE,CALCULATION.
C    Each array contains N_ARR values (data points, moments, ...)
C
        IF(I_SA.EQ.0) THEN
C
          N_ARR=N_CMP
          DO J=1,N_ARR
            X(J)=X2_CAL(J)
            IF(JFILE.EQ.1) EXPE(J)=I2_EXP(J)
            CALCULATION(J,JFILE)=I2_CAL(J)
          ENDDO
C    
        ELSEIF(I_SA.EQ.1) THEN
C
          N_ARR=N_MOM+1
          DO J=1,N_ARR
            X(J)=FLOAT(J)
            IF(JFILE.EQ.1) EXPE(J)=M_EXP(J-1)
            CALCULATION(J,JFILE)=M_CAL(J-1)
          ENDDO
C    
        ELSEIF(I_SA.EQ.2) THEN
C
          IF(METHOD.NE.'HIS') THEN
            N_ARR=N_CMP
          ELSE
            N_ARR=N_BINC
          ENDIF
          DO J=1,N_ARR
            X(J)=FLOAT(J)
            IF(JFILE.EQ.1) EXPE(J)=CH_EXP(J)
            CALCULATION(J,JFILE)=CH_CAL(J)
          ENDDO
C    
        ELSEIF(I_SA.EQ.3) THEN
C
          N_ARR=N_CMP
          DO J=1,N_ARR
            X(J)=FLOAT(J)
            IF(JFILE.EQ.1) EXPE(J)=CC_EXP(J)
            CALCULATION(J,JFILE)=CC_CAL(J)
          ENDDO
C    
        ELSEIF(I_SA.EQ.4) THEN
C
          IF(SH_AN.NE.'FOUR') THEN
            N_ARR=N_CMP
          ELSE
            N_ARR=NBIN
          ENDIF
          DO J=1,N_ARR
            X(J)=FLOAT(J)
            IF(JFILE.EQ.1) EXPE(J)=M_EXP(J)
            CALCULATION(J,JFILE)=M_CAL(J)
          ENDDO
C
        ENDIF
C
        CLOSE(NUNIT4)
        CLOSE(NUNIT5)
        CLOSE(NUNIT6)
C
      ENDDO
C
C  Computing the shift applied to all curves so that their 
C    y coordinate is always positive
C
      IF(I_POSI.EQ.1) THEN
        CALL COMPUTE_SHIFT(EXPE,CALCULATION,N_CMP,NFILE,IUO1,SHIFT)
      ELSE
        WRITE(IUO1,313)
      ENDIF
C
      NUNIT7=NU_LAST+7
      NUNIT8=NU_LAST+8
C
      DO JFILE=1,NFILE
C
C  Putting calculation file to CALC
C
        DO J=1,N_ARR
          CALC(J)=CALCULATION(J,JFILE)
        ENDDO
C
C  Shifting experiment and calculation whenever necessary
C    to ensure all similarity indices, distances, are
C              always defined
C
        IF(I_POSI.EQ.1) THEN
          IF(SHIFT.GT.SMALL) THEN
            IF(JFILE.EQ.1) CALL SHIFT_CURVE(EXPE,N_CMP,SHIFT)
            CALL SHIFT_CURVE(CALC,N_CMP,SHIFT)
          ENDIF
        ENDIF
C
        I_SW=2
        CALL NAME_FILES(I_SW,I_DIST,I_SYM,I_CUR,I_SA,NU_LAST,SH_AN,
     1                  OUTFILE(JFILE),CHFILE1,CHFILE2,CHFILE3,CHFILE4,
     2                  CHFILE5,CHFILE6,CHFILE7,CHFILE8,CHFILE9,
     3                  CHFILE10,CHFILE11)
C
C
C..........    Computation of the R-factors,         .........
C..........    similarity indices, distances,        .........
C..........    goodness of fit, kernel distances:    .........
C..........    selection of HOW it is compared       .........
C
C
C
C  The experimental and calculation arrays compared are then :
C
C          X, EXPE, CALC  varying from 1 to N_ARR 
C
C
C  Renormalizing the calculations to the experiment
C
        CALL NORMALIZE_COEF(X,EXPE,CALC,N_ARR,NFILE,I_NORM,IUO1,
     1                      OUTFILE,CNORM)
        IF(I_PRINT.EQ.1) THEN
          WRITE(IUO1,510) JFILE,CNORM
        ENDIF
C
C  Writing the renormalized experimental and calculation files
C
        DO J=1,N_ARR
C
          IF(JFILE.EQ.1) THEN
            WRITE(NUNIT7,*) X(J),EXPE(J)
          ENDIF
          WRITE(NUNIT4,*) X(J),CNORM*CALC(J)
C
        ENDDO
C
        CLOSE(NUNIT4)
        IF(JFILE.EQ.1) CLOSE(NUNIT7)
C
C  Computing the weights and writing them
C
        CALL WEIGHTS(X,EXPE,CALC,N_ARR,CNORM,IUO1,W)
        DO J=1,N_ARR
          WRITE(NUNIT8,*) X(J),W(J)
        ENDDO
C
        CLOSE(NUNIT8)
C          
        IF(I_DIST.EQ.0) THEN
C
C  Standard R-factors
C
          I_DEV=0
          CALL R_FACTOR_1(EXPE,CALC,W,N_ARR,1,CNORM,RF1,RF2,RF3,RF4,
     1                    RF5,I_DEV)
          IF(N_PAR.EQ.1) THEN
            WRITE(11,*) PAR1(JFILE),RF1
            WRITE(12,*) PAR1(JFILE),RF2
            WRITE(13,*) PAR1(JFILE),RF3
            WRITE(14,*) PAR1(JFILE),RF4
            WRITE(22,*) PAR1(JFILE),RF5
          ELSEIF(N_PAR.EQ.2) THEN
            WRITE(11,*) PAR1(JFILE),PAR2(JFILE),RF1
            WRITE(12,*) PAR1(JFILE),PAR2(JFILE),RF2
            WRITE(13,*) PAR1(JFILE),PAR2(JFILE),RF3
            WRITE(14,*) PAR1(JFILE),PAR2(JFILE),RF4
            WRITE(22,*) PAR1(JFILE),PAR2(JFILE),RF5
          ENDIF
C
C  R-factors with first order derivatives
C
          I_DEV=1
          I_FLAG=1
          CALL DERIV(EXPE,N_ARR,I_EXP_1,I_EXP_2,I_EXP_3,I_EXP_4,
     1               I_EXP_5,N_DERIV,STEP_CAL,I_FLAG)
          CALL DERIV(CALC,N_ARR,I_CAL_1,I_CAL_2,I_CAL_3,I_CAL_4,
     1               I_CAL_5,N_DERIV,STEP_CAL,I_FLAG)
          CALL R_FACTOR_1(I_EXP_1,I_CAL_1,W,N_ARR,1,CNORM,RF1,RF2,RF3,
     1                    RF4,RF5,I_DEV)
          IF(N_PAR.EQ.1) THEN
            WRITE(15,*) PAR1(JFILE),RF3
            WRITE(16,*) PAR1(JFILE),RF4
          ELSEIF(N_PAR.EQ.2) THEN
            WRITE(15,*) PAR1(JFILE),PAR2(JFILE),RF3
            WRITE(16,*) PAR1(JFILE),PAR2(JFILE),RF4
          ENDIF
C
C  Other R-factors used in medium energy ion scattering (MEIS) : 
C  D. P. Woodruff et al, Nucl. Instr. and Meth. in Phys. Res. B 183, 128 (2001)
C
          I_DEV=2
          CALL R_FACTOR_1(EXPE,CALC,W,N_ARR,1,CNORM,RF1,RF2,RF3,RF4,RF5,
     1                    I_DEV)
          IF(N_PAR.EQ.1) THEN
            WRITE(17,*) PAR1(JFILE),RF1
            WRITE(18,*) PAR1(JFILE),RF2
            WRITE(19,*) PAR1(JFILE),RF3
          ELSEIF(N_PAR.EQ.2) THEN
            WRITE(17,*) PAR1(JFILE),PAR2(JFILE),RF1
            WRITE(18,*) PAR1(JFILE),PAR2(JFILE),RF2
            WRITE(19,*) PAR1(JFILE),PAR2(JFILE),RF3
          ENDIF
C
C  Zanazzi-Jona and Pendry's R-factors
C
          I_FLAG=2
          CALL DERIV(EXPE,N_ARR,I_EXP_1,I_EXP_2,I_EXP_3,I_EXP_4,
     1               I_EXP_5,N_DERIV,STEP_CAL,I_FLAG)
          CALL DERIV(CALC,N_ARR,I_CAL_1,I_CAL_2,I_CAL_3,I_CAL_4,
     1               I_CAL_5,N_DERIV,STEP_CAL,I_FLAG)
          CALL R_FACTOR_2(EXPE,CALC,I_EXP_1,I_CAL_1,I_EXP_2,
     1                    I_CAL_2,W,N_ARR,1,CNORM,V_I,RF1,RF2)
          IF(N_PAR.EQ.1) THEN
            WRITE(20,*) PAR1(JFILE),RF1
            WRITE(21,*) PAR1(JFILE),RF2
          ELSEIF(N_PAR.EQ.2) THEN
            WRITE(20,*) PAR1(JFILE),PAR2(JFILE),RF1
            WRITE(21,*) PAR1(JFILE),PAR2(JFILE),RF2
          ENDIF
C
        ELSEIF(I_DIST.EQ.1) THEN
C
C  Similarity indices
C
          CALL SIM_INDEX(EXPE,CALC,W,N_ARR,1,CNORM,SI1,SI2,SI3,SI4,
     1                   SI5,SI6,SI7,SI8,SI9,SI10,SI11,SI12,ALPHAS,
     2                   BETAS,N_BINS)
C
          IF(N_PAR.EQ.1) THEN
            WRITE(11,*) PAR1(JFILE),SI1
            WRITE(12,*) PAR1(JFILE),SI2
            WRITE(13,*) PAR1(JFILE),SI3
            WRITE(14,*) PAR1(JFILE),SI4
            WRITE(15,*) PAR1(JFILE),SI5
            WRITE(16,*) PAR1(JFILE),SI6
            WRITE(17,*) PAR1(JFILE),SI7
            WRITE(18,*) PAR1(JFILE),SI8
            WRITE(19,*) PAR1(JFILE),SI9
            WRITE(20,*) PAR1(JFILE),SI10
            WRITE(21,*) PAR1(JFILE),SI11
            WRITE(22,*) PAR1(JFILE),SI12
          ELSEIF(N_PAR.EQ.2) THEN
            WRITE(11,*) PAR1(JFILE),PAR2(JFILE),SI1
            WRITE(12,*) PAR1(JFILE),PAR2(JFILE),SI2
            WRITE(13,*) PAR1(JFILE),PAR2(JFILE),SI3
            WRITE(14,*) PAR1(JFILE),PAR2(JFILE),SI4
            WRITE(15,*) PAR1(JFILE),PAR2(JFILE),SI5
            WRITE(16,*) PAR1(JFILE),PAR2(JFILE),SI6
            WRITE(17,*) PAR1(JFILE),PAR2(JFILE),SI7
            WRITE(18,*) PAR1(JFILE),PAR2(JFILE),SI8
            WRITE(19,*) PAR1(JFILE),PAR2(JFILE),SI9
            WRITE(20,*) PAR1(JFILE),PAR2(JFILE),SI10
            WRITE(21,*) PAR1(JFILE),PAR2(JFILE),SI11
            WRITE(22,*) PAR1(JFILE),PAR2(JFILE),SI12
          ENDIF
C
        ELSEIF(I_DIST.EQ.2) THEN
C
C  Distances
C
          CALL DISTANCE(X,EXPE,CALC,W,N_ARR,1,CNORM,DI1,DI2,DI3,
     1                  DI4,DI5,DI6,DI7,DI8,DI9,DI10,DI11,DI12,DI13,
     2                  DI14,DI15,DI16,DI17,DI18,DI19,DI20,DI21,DI22,
     3                  DI23,DI24,ALPHAD,I_BETA,SIGMAD,L)
C
          IF(N_PAR.EQ.1) THEN
            WRITE(11,*) PAR1(JFILE),DI1
            WRITE(12,*) PAR1(JFILE),DI2
            WRITE(13,*) PAR1(JFILE),DI3
            WRITE(14,*) PAR1(JFILE),DI4
            WRITE(15,*) PAR1(JFILE),DI5
            WRITE(16,*) PAR1(JFILE),DI6
            WRITE(17,*) PAR1(JFILE),DI7
            WRITE(18,*) PAR1(JFILE),DI8
            WRITE(19,*) PAR1(JFILE),DI9
            WRITE(20,*) PAR1(JFILE),DI10
            WRITE(21,*) PAR1(JFILE),DI11
            WRITE(22,*) PAR1(JFILE),DI12
            WRITE(23,*) PAR1(JFILE),DI13
            WRITE(24,*) PAR1(JFILE),DI14
            WRITE(25,*) PAR1(JFILE),DI15
            WRITE(26,*) PAR1(JFILE),DI16
            WRITE(27,*) PAR1(JFILE),DI17
            WRITE(28,*) PAR1(JFILE),DI18
            WRITE(29,*) PAR1(JFILE),DI19
            WRITE(30,*) PAR1(JFILE),DI20
            WRITE(31,*) PAR1(JFILE),DI21
            WRITE(32,*) PAR1(JFILE),DI22
            WRITE(33,*) PAR1(JFILE),DI23
            WRITE(34,*) PAR1(JFILE),DI24
          ELSEIF(N_PAR.EQ.2) THEN
            WRITE(11,*) PAR1(JFILE),PAR2(JFILE),DI1
            WRITE(12,*) PAR1(JFILE),PAR2(JFILE),DI2
            WRITE(13,*) PAR1(JFILE),PAR2(JFILE),DI3
            WRITE(14,*) PAR1(JFILE),PAR2(JFILE),DI4
            WRITE(15,*) PAR1(JFILE),PAR2(JFILE),DI5
            WRITE(16,*) PAR1(JFILE),PAR2(JFILE),DI6
            WRITE(17,*) PAR1(JFILE),PAR2(JFILE),DI7
            WRITE(18,*) PAR1(JFILE),PAR2(JFILE),DI8
            WRITE(19,*) PAR1(JFILE),PAR2(JFILE),DI9
            WRITE(20,*) PAR1(JFILE),PAR2(JFILE),DI10
            WRITE(21,*) PAR1(JFILE),PAR2(JFILE),DI11
            WRITE(22,*) PAR1(JFILE),PAR2(JFILE),DI12
            WRITE(23,*) PAR1(JFILE),PAR2(JFILE),DI13
            WRITE(24,*) PAR1(JFILE),PAR2(JFILE),DI14
            WRITE(25,*) PAR1(JFILE),PAR2(JFILE),DI15
            WRITE(26,*) PAR1(JFILE),PAR2(JFILE),DI16
            WRITE(27,*) PAR1(JFILE),PAR2(JFILE),DI17
            WRITE(28,*) PAR1(JFILE),PAR2(JFILE),DI18
            WRITE(29,*) PAR1(JFILE),PAR2(JFILE),DI19
            WRITE(30,*) PAR1(JFILE),PAR2(JFILE),DI20
            WRITE(31,*) PAR1(JFILE),PAR2(JFILE),DI21
            WRITE(32,*) PAR1(JFILE),PAR2(JFILE),DI22
            WRITE(33,*) PAR1(JFILE),PAR2(JFILE),DI23
            WRITE(34,*) PAR1(JFILE),PAR2(JFILE),DI24
          ENDIF
C
        ELSEIF(I_DIST.EQ.3) THEN
C
C  Goodness of fit
C
          CALL GOODNESS(EXPE,CALC,W,N_ARR,1,CNORM,GF1,GF2,GF3,
     1                  GF4,GF5,GF6,GF7,GF8,GF9,GF10,GF11,GF12,
     2                  N_BING,ALPHAG)
C
          IF(N_PAR.EQ.1) THEN
            WRITE(11,*) PAR1(JFILE),GF1
            WRITE(12,*) PAR1(JFILE),GF2
            WRITE(13,*) PAR1(JFILE),GF3
            WRITE(14,*) PAR1(JFILE),GF4
            WRITE(15,*) PAR1(JFILE),GF5
            WRITE(16,*) PAR1(JFILE),GF6
            WRITE(17,*) PAR1(JFILE),GF7
            WRITE(18,*) PAR1(JFILE),GF8
            WRITE(19,*) PAR1(JFILE),GF9
            WRITE(20,*) PAR1(JFILE),GF10
            WRITE(21,*) PAR1(JFILE),GF11
            WRITE(22,*) PAR1(JFILE),GF12
          ELSEIF(N_PAR.EQ.2) THEN
            WRITE(11,*) PAR1(JFILE),PAR2(JFILE),GF1
            WRITE(12,*) PAR1(JFILE),PAR2(JFILE),GF2
            WRITE(13,*) PAR1(JFILE),PAR2(JFILE),GF3
            WRITE(14,*) PAR1(JFILE),PAR2(JFILE),GF4
            WRITE(15,*) PAR1(JFILE),PAR2(JFILE),GF5
            WRITE(16,*) PAR1(JFILE),PAR2(JFILE),GF6
            WRITE(17,*) PAR1(JFILE),PAR2(JFILE),GF7
            WRITE(18,*) PAR1(JFILE),PAR2(JFILE),GF8
            WRITE(19,*) PAR1(JFILE),PAR2(JFILE),GF9
            WRITE(20,*) PAR1(JFILE),PAR2(JFILE),GF10
            WRITE(21,*) PAR1(JFILE),PAR2(JFILE),GF11
            WRITE(22,*) PAR1(JFILE),PAR2(JFILE),GF12
          ENDIF
C
        ELSEIF(I_DIST.EQ.4) THEN
C
C  Kernel distances
C
          IF(ABS(SIGMAD).LT.SMALL) THEN
            SIGMAD= 1.0
            WRITE(IUO1,311)
          ENDIF
C
          CALL KERNEL(EXPE,CALC,W,N_ARR,ALPHAK,SIGMAK,L,CNORM,KD1,KD2,
     1                KD3,KD4,KD5,KD6,KD7,KD8,KD9,KD10,KD11,KD12)
C
          IF(N_PAR.EQ.1) THEN
            WRITE(11,*) PAR1(JFILE),KD1
            WRITE(12,*) PAR1(JFILE),KD2
            WRITE(13,*) PAR1(JFILE),KD3
            WRITE(14,*) PAR1(JFILE),KD4
            WRITE(15,*) PAR1(JFILE),KD5
            WRITE(16,*) PAR1(JFILE),KD6
            WRITE(17,*) PAR1(JFILE),KD7
            WRITE(18,*) PAR1(JFILE),KD8
            WRITE(19,*) PAR1(JFILE),KD9
            WRITE(20,*) PAR1(JFILE),KD10
            WRITE(21,*) PAR1(JFILE),KD11
            WRITE(22,*) PAR1(JFILE),KD12
          ELSEIF(N_PAR.EQ.2) THEN
            WRITE(11,*) PAR1(JFILE),PAR2(JFILE),KD1
            WRITE(12,*) PAR1(JFILE),PAR2(JFILE),KD2
            WRITE(13,*) PAR1(JFILE),PAR2(JFILE),KD3
            WRITE(14,*) PAR1(JFILE),PAR2(JFILE),KD4
            WRITE(15,*) PAR1(JFILE),PAR2(JFILE),KD5
            WRITE(16,*) PAR1(JFILE),PAR2(JFILE),KD6
            WRITE(17,*) PAR1(JFILE),PAR2(JFILE),KD7
            WRITE(18,*) PAR1(JFILE),PAR2(JFILE),KD8
            WRITE(19,*) PAR1(JFILE),PAR2(JFILE),KD9
            WRITE(20,*) PAR1(JFILE),PAR2(JFILE),KD10
            WRITE(21,*) PAR1(JFILE),PAR2(JFILE),KD11
            WRITE(22,*) PAR1(JFILE),PAR2(JFILE),KD12
          ENDIF
C
        ENDIF
C
        CLOSE(NUNIT3)
C
      ENDDO
C
      DO JR=1,N_RF
        NUNIT=JR+10
        CLOSE(NUNIT)
      ENDDO
C
C
      IF(I_DIST.EQ.0) THEN 
        WRITE(IUO1,103) 
      ELSEIF(I_DIST.EQ.1) THEN
        WRITE(IUO1,104) 
      ELSEIF(I_DIST.EQ.2) THEN
        WRITE(IUO1,105) 
      ELSEIF(I_DIST.EQ.3) THEN
        WRITE(IUO1,106) 
      ELSEIF(I_DIST.EQ.4) THEN
        WRITE(IUO1,107) 
      ENDIF
C
CST      CLOSE(IUO1)
C
   1  FORMAT(A40)
   2  FORMAT(10X,17A4) 
   3  FORMAT(9X,A1)  
  11  FORMAT(//,10X,'<<<<<  N_FILES SHOULD BE LARGER THAN ',I5,
     1      '  >>>>>',//)
  14  FORMAT(8X,I2,9X,I1,9X,I1,9X,I1)
  15  FORMAT(9X,I1,6X,F7.2,6X,I1)
  16  FORMAT(9X,I1,9X,I1,9X,I1,9X,I1)
  17  FORMAT(9X,I1,6X,F7.2,3X,F7.2,3X,F7.2)
  18  FORMAT(9X,I1,4X,F9.2)
  19  FORMAT(8X,F6.3)
  20  FORMAT(8X,F6.3,4X,F6.3,2X,I4)
  21  FORMAT(8X,F6.3,5X,I1,8X,I2,6X,F7.2)
  22  FORMAT(8X,F6.3,4X,I2,6X,F7.2)
  23  FORMAT(5X,I5,5X,I5,6X,A4,9X,I1)
  24  FORMAT(9X,I1,6X,D7.2,3X,D7.2)
  25  FORMAT(9X,I1,7X,A3,8X,I2,7X,I3)
  26  FORMAT(9X,I1,6X,F7.2)
  27  FORMAT(5X,I5,5X,I5,6X,A4,9X,I1)
  28  FORMAT(9X,I1)
  29  FORMAT(6X,I4,6X,F7.2)
  34  FORMAT(9X,A40)
  35  FORMAT(9X,A40,1X,F13.6,1X,F13.6)
  36  FORMAT(6X,'|',/,6X,'CONTOUR SHAPE DESCRIPTOR WRITTEN IN ',A48,//)
  47  FORMAT(6X,'CHAIN CODE WRITTEN IN ',A48,//)
  48  FORMAT(6X,'CHORDS WRITTEN IN ',A48,//)
  49  FORMAT(6X,'MOMENTS WRITTEN IN ',A48)
  51  FORMAT(39X,A1)
C
 100  FORMAT(/////,'******************************', 
     1       '****************************************************')  
 101  FORMAT('********',66X,'********')     
 102  FORMAT('********  ',15A4,'    ********') 
 103  FORMAT(////,'*********************   END OF R-FACTOR ANALYSIS   ',
     1       '*********************',/////) 
   
 104  FORMAT(////,'*********************   END OF SIMILARITY INDEX ',
     1       'ANALYSIS   *********************',/////) 
 105  FORMAT(////,'*********************   END OF MATHEMATICAL DISTAN',
     1       'CE ANALYSIS   *********************',/////) 
 106  FORMAT(////,'*********************   END OF GOODNESS OF FIT ',
     1       'ANALYSIS   *********************',/////) 
 107  FORMAT(////,'*********************   END OF KERNEL DISTANCE ',
     1       'ANALYSIS   *********************',/////) 
 114  FORMAT(16X,I2,9X,I1,9X,I1,9X,I1,9X,'N_PAR,NORM,I_SCALE,I_NORM')       
 115  FORMAT(17X,I1,6X,F7.2,6X,I1,19X,'I_SYM,SYM,I_POSI')   
 116  FORMAT(17X,I1,9X,I1,9X,I1,9X,I1,9X,'I_DIST,I_CUR,I_SA,I_PRINT',
     1        ///)
 117  FORMAT(17X,I1,6X,F7.2,3X,F7.2,3X,F7.2,6X,'I_WEIGHT,ALPHA,BETA,',
     1       'SIGMA')
 118  FORMAT(17X,I1,4X,F9.2,26X,'I_SHIFT,MAXW',///)
 119  FORMAT(16X,F6.3,35X,'V_I',///)
 120  FORMAT(16X,F6.3,4X,F6.3,2X,I4,19X,'ALPHAS,BETAS,N_BINS',///)
 121  FORMAT(16X,F6.3,5X,I1,8X,I2,6X,F7.2,6X,'ALPHAD,I_BETA,L,SIGMAD',
     1       ///)
 122  FORMAT(16X,F6.3,4X,I2,6X,F7.2,16X,'ALPHAK,L,SIGMAK',///)
 123  FORMAT(13X,I5,5X,I5,6X,A4,19X,'N_GRID,N_MOM,BASIS')
 124  FORMAT(17X,I1,6X,F7.2,3X,F7.2,16X,'I_ALG,MU,NU',///)
 125  FORMAT(17X,I1,7X,A3,8X,I2,7X,I3,9X,'I_CHORD,METHOD,VALUE,N_BINC',
     1       /)
 126  FORMAT(17X,I1,6X,F7.2,26X,'N_CONNECT,SCALEC',///)
 127  FORMAT(13X,I5,5X,I5,6X,A4,9X,I1,9X,'N_BIN,N_LEN,SH_AN,I_FOU',///)
 128  FORMAT(17X,I1,39X,'INORM',///)
 129  FORMAT(14X,I4,6X,F7.2,26X,'N_BING,ALPHAG',///)
 140  FORMAT(6X,'CALCULATION PARAMETERS : GENERAL',/)
 141  FORMAT(6X,'CALCULATION PARAMETERS : WEIGHTS',/)
 142  FORMAT(6X,'CALCULATION PARAMETERS : R-FACTORS',/)
 143  FORMAT(6X,'CALCULATION PARAMETERS : SIMILARITY INDICES',/)
 144  FORMAT(6X,'CALCULATION PARAMETERS : DISTANCES',/)
 145  FORMAT(6X,'CALCULATION PARAMETERS : KERNEL DISTANCES',/)
 146  FORMAT(6X,'CALCULATION PARAMETERS : MOMENTS',/)
 147  FORMAT(6X,'CALCULATION PARAMETERS : CHORDS',/)
 148  FORMAT(6X,'CALCULATION PARAMETERS : CHAIN CODES',/)
 149  FORMAT(6X,'CALCULATION PARAMETERS : CONTOUR',/)
 150  FORMAT(6X,'EXPERIMENTAL FILE : ',A40,' PARAMETER 1 :',//)      
 151  FORMAT(6X,'CALCULATION FILES : ',A40,1X,F13.6)
 152  FORMAT(26X,A40,1X,F13.6,1X,F13.6)                                         
 153  FORMAT(6X,'EXPERIMENTAL FILE : ',A40,' PARAMETER 1 :', 1X,
     1       'PARAMETER 2 :', //)      
 154  FORMAT(6X,'CALCULATION FILES : ',A40,1X,F13.6,1X,F13.6)
 155  FORMAT(6X,'CALCULATION PARAMETERS : GOODNESS OF FIT',/)
C
 203  FORMAT('**************************************************',      
     1       '********************************',//////////)  
C
 308  FORMAT(//,10X,'-->  DIMENSIONING ERROR : NMAX SHOULD BE ',
     1      'AT LEAST ',I9,//)
 309  FORMAT(//,10X,'-->  INPUT DATA ERROR : N_LEN SHOULD BE ',
     1      'LESS THAN ',I9,//)
 310  FORMAT(//,10X,'-->  INPUT DATA ERROR : N_BIN SHOULD BE ',
     1      'LESS THAN ',I9,//)
 311  FORMAT(//,10X,'-->  TO AVOID DIVERGENCE IN THE KERNEL FUNCTIONS,',
     1       'C HAS BEEN SET TO 1.0',//)
 312  FORMAT(//,10X,'----> STEP USED FOR THE OVERALL COMPARISON : ',
     1       F8.2,//)
 313  FORMAT(//,10X,'-->  NO AUTOMATIC SHIFTING OF Y COORDINATES: ',
     1          'THIS MIGHT LEAD ',/,10X,
     2       '-->  TO UNDEFINED COMPARISONS FOR SOME METHODS !',//)
C
 400  FORMAT(//,6X,'----> R-FACTOR ANALYSIS:')
 401  FORMAT(//,6X,'----> SIMILARITY INDEX ANALYSIS:')
 402  FORMAT(//,6X,'----> MATHEMATICAL DISTANCE ANALYSIS:')
 403  FORMAT(//,6X,'----> GOODNESS OF FIT ANALYSIS:')
 404  FORMAT(//,6X,'----> KERNEL DISTANCE ANALYSIS:')
 405  FORMAT(20X,'DESCRIPTOR USED: MODULATION FUNCTION (CHI)',//)
 406  FORMAT(20X,'DESCRIPTOR USED: COORDINATES OF CURVE',//)
 407  FORMAT(20X,'DESCRIPTOR USED: CONTINUOUS GEOMETRIC MOMENTS',//)
 408  FORMAT(20X,'DESCRIPTOR USED: CONTINUOUS LEGENDRE MOMENTS',//)
 409  FORMAT(20X,'DESCRIPTOR USED: DISCRETE CHEBYSHEV MOMENTS',//)
 410  FORMAT(20X,'DESCRIPTOR USED: DISCRETE KRAWTCHOUK MOMENTS',//)
 411  FORMAT(20X,'DESCRIPTOR USED: DISCRETE HAHN MOMENTS',//)
 412  FORMAT(20X,'DESCRIPTOR USED: DISCRETE MEIXNER MOMENTS',//)
 413  FORMAT(20X,'DESCRIPTOR USED: DISCRETE CHARLIER MOMENTS',//)
 414  FORMAT(20X,'DESCRIPTOR USED: DISCRETE SHMALIY MOMENTS',//)
 415  FORMAT(20X,'CHORD DESCRIPTOR USED: CHORD LENGTH FROM POINT',//)
 416  FORMAT(20X,'CHORD DESCRIPTOR USED: ARC-CHORD DISTANCE',//)
 417  FORMAT(20X,'CHORD DESCRIPTOR USED: CHORD LENGTH ALONG ANGLE',//)
 418  FORMAT(20X,'CHAIN CODE DESCRIPTOR CONNECTIVITY: 3',//)
 419  FORMAT(20X,'CHAIN CODE DESCRIPTOR CONNECTIVITY: 5',//)
 420  FORMAT(20X,'CHAIN CODE DESCRIPTOR CONNECTIVITY: 9',//)
 421  FORMAT(20X,'DESCRIPTOR USED: CUMULATIVE DISTRIBUTION FUNCTION',//)
 422  FORMAT(20X,'DESCRIPTOR USED: CURVATURE OF CURVE',//)
 423  FORMAT(20X,'DESCRIPTOR USED: CUI-FEMIANI-HU-WONKA-RAZDAN ',
     1           'SIGNATURE OF CURVE',//)
 437  FORMAT(20X,'SHAPE DESCRIPTOR USED: CENTROID DISTANCE',//)
 438  FORMAT(20X,'SHAPE DESCRIPTOR USED: TANGENT ANGLE',//)
 439  FORMAT(20X,'SHAPE DESCRIPTOR USED: CURVATURE FUNCTION',//)
 440  FORMAT(20X,'SHAPE DESCRIPTOR USED: TRIANGLE AREA',//)
 441  FORMAT(20X,'SHAPE DESCRIPTOR USED: BEAM ANGLE STATISTICS',//)
 442  FORMAT(20X,'SHAPE DESCRIPTOR USED: 8-CONNECTIVITY CHAIN CODE',//)
 443  FORMAT(20X,'SHAPE DESCRIPTOR USED: CHORD LENGTH',//)
 444  FORMAT(20X,'SHAPE DESCRIPTOR USED: CHORD ANGLE',//)
 445  FORMAT(20X,'SHAPE DESCRIPTOR USED: ARCH CHORD DISTANCE',//)
 446  FORMAT(20X,'SHAPE DESCRIPTOR USED: FOURIER COEFFICIENTS',//)
 450  FORMAT(//,6X,'----> PRIOR NORMALIZATION OF CURVES: NONE')
 451  FORMAT(//,6X,'----> PRIOR NORMALIZATION OF CURVES: SECOND ORDER ',
     1       'CENTRAL MOMENT = 1')
 452  FORMAT(//,6X,'----> PRIOR NORMALIZATION OF CURVES: TO UNIT AREA')
 453  FORMAT(//,6X,'----> PRIOR NORMALIZATION OF CURVES: TO MAXIMUM ')
 454  FORMAT(//,6X,'----> PRIOR NORMALIZATION OF CURVES: DECIMAL ',
     1       'SCALING')
 455  FORMAT(//,6X,'----> PRIOR NORMALIZATION OF CURVES: TO [0,1]')
 456  FORMAT(//,6X,'----> SCALING THE EACH CALCULATION TO MIN-MAX OF ',
     1       'EXPERIMENT PRIOR TO ANALYSIS',//)
 457  FORMAT(6X,'----> SHAPE ANALYSIS: I_NORM HAS BEEN RESET TO ',
     1       'ZERO',//)
 458  FORMAT(//,6X,'----> GOODNESS OF FIT ANALYSIS: I_POSI ',
     1             'AUTOMATICALLY SET TO 1')
 459  FORMAT(20X,'----> CUMULATIVE DISTRIBUTION FUNCTION: NORM ',
     1             'AUTOMATICALLY SET TO 0',//)
 460  FORMAT(20X,'----> CUMULATIVE DISTRIBUTION FUNCTION: I_SCALE ',
     1             'AUTOMATICALLY SET TO 1',//)
 461  FORMAT(20X,'----> CONTOUR-BASED DESCRIPTOR: NORM AUTOMATICALLY ',
     1           'SET TO 0 BECAUSE',/,52X,'CONTOUR WILL BE ITSELF ',
     2           'NORMALIZED',//)
 462  FORMAT(//,10X,'----> I_DIST PARAMETER EXTERNALLY SET TO ',I1)
 463  FORMAT(//,10X,'----> I_SA PARAMETER EXTERNALLY SET TO ',I1)
 464  FORMAT(//,10X,'----> EXTERNAL INFORMATION READ FROM THE SCRIPT ',
     1              'FILE proc_nrfactor')
C
 510  FORMAT(10X,'----> NORMALIZATION COEFFICIENT FOR CALCULATION ',I4,
     1       ':  ',E12.6)
C
      END
C
C=======================================================================
C
      SUBROUTINE CHECK_CALC_FILE(CALCFILE,NFILE,IUO1,STEP_CAL)
C
C  This subroutine reads each calculated file to compute its lower 
C     bound, upper bound and step. They are stored. Then, they are 
C     all compared to ensure that they are equal for each file.
C     Otherwise, the code is stopped with an error message.
C
C  Input parameters:
C
C       CALCFILE :  name of the calculation files
C       NFILE    :  number of calculation files
C       IUO1     :  checkfile index for printing
C
C
C  Input parameters:
C
C       STEP_CAL :  step for the calculation files
C
C   Author :  D. Sébilleau
C  
C                                           Last modified : 9 Sep 2014
C
      PARAMETER (N_SIZE=1000,N_FILES=1000)
C
      REAL*4 X_CAL(N_SIZE),I_CAL,X,Y
      REAL*4 DIFF,STEP1,STEP2
      REAL*4 LB_MIN,LB_MAX,UB_MIN,UB_MAX,STEP_MIN,STEP_MAX
      REAL*4 L_BOUND(N_FILES),U_BOUND(N_FILES),STEP(N_FILES)
C
      INTEGER LINES(N_SIZE),LINE_MIN,LINE_MAX,FILE_STOP
C
      CHARACTER*40 CALCFILE(N_FILES)
C
      DATA COMP /0.1/
C
      FILE_STOP=0
      IRET=0
C
C  Checking if dimensioning N_SIZE is large enough
C
      DO JFILE=1,NFILE
C
        OPEN(UNIT=1, FILE=CALCFILE(JFILE), STATUS='unknown')
C
        N_PO=0
        DO JLINE=1,N_SIZE+5
          READ(1,*,END=25) X,Y
          N_PO=N_PO+1
        ENDDO
  25    IF(N_PO.GT.N_SIZE) THEN 
          IRET=1
          GOTO 15
        ENDIF
C
        CLOSE(1)
C
      ENDDO
C
C  Computing the lower bound, upper bound, number of lines  
C    and step for each calculation files
C
      DO JFILE=1,NFILE
C
        OPEN(UNIT=1, FILE=CALCFILE(JFILE), STATUS='unknown')
C
        STEP1=9999999.0
        STEP2=0.
        N_PO=0
        DO JLINE=1,N_SIZE
          READ(1,*,END=5) X_CAL(JLINE),I_CAL
          N_PO=N_PO+1
          IF(JLINE.GE.2) THEN
            DIFF=X_CAL(JLINE)-X_CAL(JLINE-1)
            STEP1=MIN(DIFF,STEP1)
            STEP2=MAX(DIFF,STEP2)
          ENDIF
        ENDDO
C
   5    CONTINUE
        LINES(JFILE)=N_PO
        L_BOUND(JFILE)=X_CAL(1)
        U_BOUND(JFILE)=X_CAL(N_PO)
        IF(ABS(STEP1-STEP2).LT.COMP) THEN
          STEP(JFILE)=STEP1
        ELSE
          IRET=2
          FILE_STOP=JFILE
        ENDIF
C
        CLOSE(1)
C
      ENDDO
C
      IF(IRET.EQ.2) GOTO 15
C
      LB_MIN=9999999.0
      LB_MAX=-9999999.0
      UB_MIN=9999999.0
      UB_MAX=-9999999.0
      STEP_MIN=9999999.0
      STEP_MAX=0.
      LINE_MIN=9999999
      LINE_MAX=0
C
C  Computing the lower bound, upper bound number of lines  
C    and step for the set of calculation files
C
      DO JFILE=1,NFILE
C
        LB_MIN=MIN(L_BOUND(JFILE),LB_MIN)
        LB_MAX=MAX(L_BOUND(JFILE),LB_MAX)
        UB_MIN=MIN(U_BOUND(JFILE),UB_MIN)
        UB_MAX=MAX(U_BOUND(JFILE),UB_MAX)
        STEP_MIN=MIN(STEP(JFILE),STEP_MIN)
        STEP_MAX=MAX(STEP(JFILE),STEP_MAX)
        LINE_MIN=MIN(LINES(JFILE),LINE_MIN)
        LINE_MAX=MAX(LINES(JFILE),LINE_MAX)
C
      ENDDO
C
C  Inconsistencies between the different files
C
      IF(ABS(LB_MIN-LB_MAX).GT.COMP) THEN
        IRET=3
      ENDIF
      IF(ABS(UB_MIN-UB_MAX).GT.COMP) THEN
        IRET=4
      ENDIF
      IF(ABS(STEP_MIN-STEP_MAX).GT.COMP) THEN
        IRET=5
      ENDIF
      IF(ABS(LINE_MIN-LINE_MAX).GT.0) THEN
        IRET=6
      ENDIF
C
C  Calculation step
C
      IF(IRET.EQ.0) THEN
        STEP_CAL=STEP_MIN
      ENDIF
C
C  Stops for inconsistencies
C
  15  IF(IRET.EQ.1) THEN
        WRITE(IUO1,10) 
        STOP
      ELSEIF(IRET.EQ.2) THEN
        WRITE(IUO1,20) FILE_STOP
        STOP
      ELSEIF(IRET.EQ.3) THEN
        WRITE(IUO1,30) 
        STOP
      ELSEIF(IRET.EQ.4) THEN
        WRITE(IUO1,40) 
        STOP
      ELSEIF(IRET.EQ.5) THEN
        WRITE(IUO1,50) 
        STOP
      ELSEIF(IRET.EQ.6) THEN
        WRITE(IUO1,60) 
        STOP
      ENDIF
C
C  Format
C
  10  FORMAT(//,10X,'<<<<<  N_SIZE NOT SUFFICIENT ! >>>>>',/,
     1          10X,'<<<<<  INCREASE IT EVERYWHERE  >>>>>',//)
  20  FORMAT(//,10X,'<<<<<  IRREGULAR STEP IN FILE ',I4,'  >>>>>',//)
  30  FORMAT(//,10X,'<<<<<  LOWER BOUNDS DIFFERENT FOR X  >>>>>',/,
     1          10X,'<<<<<   IN AT LEAST TWO CALC FILES   >>>>>',//)
  40  FORMAT(//,10X,'<<<<<  UPPER BOUNDS DIFFERENT FOR X  >>>>>',/,
     1          10X,'<<<<<   IN AT LEAST TWO CALC FILES   >>>>>',//)
  50  FORMAT(//,10X,'<<<<<  STEP DIFFERENT FOR X   >>>>>',/,
     1          10X,'<<<<<  IN AT LEAST TWO FILES  >>>>>',//)
  60  FORMAT(//,10X,'<<<<<  NUMBER OF LINES DIFFERENT  >>>>>',/,
     1          10X,'<<<<<    IN AT LEAST TWO FILES    >>>>>',//)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE NAME_FILES(I_SW,I_DIST,I_SYM,I_CUR,I_SA,NU_LAST,SH_AN,
     1                      INFILE,CHFILE1,CHFILE2,CHFILE3,CHFILE4,
     2                      CHFILE5,CHFILE6,CHFILE7,CHFILE8,CHFILE9,
     3                      CHFILE10,CHFILE11)
C
C  This subroutine names interpolated/symmetrized experimental files,
C    and the transformed files (moments, ...)
C
C  Input parameters:
C
C       I_SW:
C           = 1 :  experimental file
C           = 2 :  calculation file   
C 
C       I_DIST  :  R-factor/ Similarity index ... switch
C       I_SYM   :  symmetrization switch (experimental file only)
C       I_CUR   :  switch for calculation of modulation function
C       I_SA    :  shape analysis switch
C       NU_LAST :  index of last file opened
C       INFILE  :  name of the input file
C
C  Input parameters:
C
C       CHFILEn :  name of the different output files
C
C
C   Author :  D. Sébilleau
C
C                                            Last modified : 16 Jan 2015
C
      CHARACTER*4  SH_AN
      CHARACTER*40 INFILE
      CHARACTER*48 CHFILE1,CHFILE2,CHFILE3,CHFILE4,CHFILE5,CHFILE6
      CHARACTER*48 CHFILE7,CHFILE9,CHFILE10,CHFILE11
      CHARACTER*52 CHFILE8
C
      NUNIT3=NU_LAST+3
      NUNIT4=NU_LAST+4
      NUNIT5=NU_LAST+5
      NUNIT6=NU_LAST+6
      NUNIT7=NU_LAST+7
      NUNIT8=NU_LAST+8
      NUNIT9=NU_LAST+9
      NUNIT10=NU_LAST+10
C
C  Reducing the name of the input file to its exact size.
C    Finding the position of slash and dot in the file name
C
      N_DOT=1
      DO J_CHAR=1,40
        IF(INFILE(J_CHAR:J_CHAR).EQ.'.') GOTO 10
        N_DOT=N_DOT+1
      ENDDO
  10  CONTINUE
C
      N_CHAR=0
      DO J_CHAR=1,40
        IF(INFILE(J_CHAR:J_CHAR).EQ.' ') GOTO 20
        N_CHAR=N_CHAR+1
      ENDDO
  20  CONTINUE
C
      N_SL=1
      DO J_CHAR=1,40
        IF(INFILE(J_CHAR:J_CHAR).EQ.'/') GOTO 30
        N_SL=N_SL+1
      ENDDO
   30 CONTINUE
C
      IF(I_SW.EQ.1) THEN
C
C  Input: experimental file case
C
        CHFILE1=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_ren'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=NUNIT7, FILE=CHFILE1, STATUS='unknown')
        CHFILE2=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_int'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=NUNIT3, FILE=CHFILE2, STATUS='unknown')
        IF(I_SYM.NE.0) THEN
        CHFILE3=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_sym'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=NUNIT4, FILE=CHFILE3, STATUS='unknown')
        ENDIF
C
        IF(I_CUR.EQ.1) THEN
        CHFILE5=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_chi'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=NUNIT6, FILE=CHFILE5, STATUS='unknown')
        ENDIF
        IF((I_CUR.EQ.2).OR.(I_DIST.EQ.3)) THEN
        CHFILE4=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_cdf'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=NUNIT5, FILE=CHFILE4, STATUS='unknown')
        ENDIF
        IF(I_CUR.EQ.3) THEN
        CHFILE10=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//
     1           '_cur'//INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=NUNIT9, FILE=CHFILE10, STATUS='unknown')
        ENDIF
        IF(I_CUR.EQ.4) THEN
        CHFILE11=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//
     1           '_cfh'//INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=NUNIT10, FILE=CHFILE11, STATUS='unknown')
        ENDIF
C
      ELSE
C
C  Input: calculation file case
C
        CHFILE2=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_ren'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=NUNIT4, FILE=CHFILE2, STATUS='unknown')
C
        IF(I_CUR.EQ.1) THEN
        CHFILE5=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_chi'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=NUNIT6, FILE=CHFILE5, STATUS='unknown')
        ENDIF
        IF((I_CUR.EQ.2).OR.(I_DIST.EQ.3)) THEN
        CHFILE4=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_cdf'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=NUNIT5, FILE=CHFILE4, STATUS='unknown')
        ENDIF
        IF(I_CUR.EQ.3) THEN
        CHFILE10=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//
     1           '_cur'//INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=NUNIT9, FILE=CHFILE10, STATUS='unknown')
        ENDIF
        IF(I_CUR.EQ.4) THEN
        CHFILE11=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//
     1           '_cfh'//INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=NUNIT10, FILE=CHFILE11, STATUS='unknown')
        ENDIF
C
        CHFILE9=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_wgt'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=NUNIT8, FILE=CHFILE9, STATUS='unknown')
C
      ENDIF
C
C  Input: all cases
C
      IF(I_SA.EQ.1) THEN
        CHFILE6=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_rec'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=99, FILE=CHFILE6, STATUS='unknown')
        CHFILE7=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_mom'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=98, FILE=CHFILE7, STATUS='unknown')
      ELSEIF(I_SA.EQ.2) THEN
        CHFILE7=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_cho'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=98, FILE=CHFILE7, STATUS='unknown')
      ELSEIF(I_SA.EQ.3) THEN
        CHFILE7=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_cco'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=98, FILE=CHFILE7, STATUS='unknown')
      ELSEIF(I_SA.EQ.4) THEN
        CHFILE6=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_ctr'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=99, FILE=CHFILE6, STATUS='unknown')
        CHFILE7=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//'_sha'//
     1          INFILE(N_DOT:N_CHAR)
        OPEN(UNIT=98, FILE=CHFILE7, STATUS='unknown')
        IF(SH_AN.EQ.'FOUR') THEN
          CHFILE8=INFILE(1:N_SL-1)//'/div'//INFILE(N_SL:N_DOT-1)//
     1            '_ctr_rec'//INFILE(N_DOT:N_CHAR)
          OPEN(UNIT=97, FILE=CHFILE8, STATUS='unknown')
        ENDIF
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE FILE_INTERSECTION(NU_LAST,EXPFILE,CALFILE,IUO1,
     1                             N_EXP,N_CAL,STEP_EXP,STEP_CAL,N_CMP,
     2                             J_IN,J_FI,X_EXP,I_EXP,X_CAL,I_CAL)
C
C  This subroutine reads the experimental and calculation files
C    and work out their intersection
C
C  Input parameters:
C
C       NU_LAST  :  last file index used before calling subroutine
C       EXPFILE  :  experimental filename
C       CALCFILE :  calculation filename
C       IUO1     :  checkfile index for printing
C
C  Output parameters:
C
C       N_EXP    :  number of points in experimental curve
C       N_CAL    :  number of points in calculation curve
C       STEP_EXP :  step of experimental file
C       STEP_CAL :  step of calculation file
C       N_CMP    :  number of points in intersection
C       J_IN     :  position of first point of intersection in 
C                   calculation grid
C       J_FI     :  position of last point of intersection in 
C                   calculation grid
C       X_EXP    :  x coordinates of experimental file
C       I_EXP    :  y coordinates of experimental file 
C       X_CAL    :  x coordinates of calculation file
C       I_CAL    :  y coordinates of calculation file 
C
C
C   Author :  D. Sébilleau
C  
C                                           Last modified :  5 Sep 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 X_CAL(N_SIZE),I_CAL(N_SIZE)
      REAL*4 X_EXP(N_SIZE),I_EXP(N_SIZE)
C
      CHARACTER*40 EXPFILE,CALFILE
C
      DATA COMP  /0.1/
C
C    Reading the experimental file and checking the step
C
      NU_EXP=NU_LAST+1
      OPEN(UNIT=NU_EXP, FILE=EXPFILE, STATUS='unknown')
C  
      N_EXP=0
      STEP1=9999999.0
      STEP2=0.
C       
      DO JLINE=1,N_SIZE
        READ(NU_EXP,*,END=15) X_EXP(JLINE),I_EXP(JLINE)
        N_EXP=N_EXP+1
        IF(JLINE.GE.2) THEN
          STEP=X_EXP(JLINE)-X_EXP(JLINE-1)
          STEP1=MIN(STEP,STEP1)
          STEP2=MAX(STEP,STEP2)
        ENDIF
      ENDDO
C
  15  WRITE(IUO1,20) EXPFILE,N_EXP,X_EXP(1),X_EXP(N_EXP)
      IF(ABS(STEP1-STEP2).LT.COMP) THEN
        STEP_EXP=STEP1
        WRITE(IUO1,30) STEP_EXP
      ELSE
        WRITE(IUO1,40)
      ENDIF
C
      CLOSE(NU_EXP)
C
      IF(N_EXP.GT.N_SIZE) THEN
        WRITE(IUO1,10) N_EXP
        STOP
      ENDIF
C
C    Reading the first calculation file and checking the step
C
      NU_CAL=NU_LAST+2
      OPEN(UNIT=NU_CAL, FILE=CALFILE, STATUS='unknown')
C
      N_CAL=0
      STEP1=9999999.0
      STEP2=0.
C       
      DO JLINE=1,N_SIZE
        READ(NU_CAL,*,END=25) X_CAL(JLINE),I_CAL(JLINE)
        N_CAL=N_CAL+1
        IF(JLINE.GE.2) THEN
          STEP=X_CAL(JLINE)-X_CAL(JLINE-1)
          STEP1=MIN(STEP,STEP1)
          STEP2=MAX(STEP,STEP2)
        ENDIF
      ENDDO
  25  IF(ABS(STEP1-STEP2).LT.COMP) THEN
        STEP_CAL=STEP1
      ELSE
        WRITE(IUO1,40)
      ENDIF
      WRITE(IUO1,20) CALFILE,N_CAL,X_CAL(1),X_CAL(N_CAL)
      WRITE(IUO1,30) STEP_CAL
C
      CLOSE(NU_CAL)
C
C..........  Intersection of the experimental    ..........
C..........      and the calculation grids       ..........
C
      X_IN=MAX(X_EXP(1),X_CAL(1))
      X_FI=MIN(X_EXP(N_EXP),X_CAL(N_CAL))
C
C    Position of these bounds in the calculation grid 
C            (used as the computing grid)
C
      CALL LOCATE(X_CAL,N_CAL,X_IN,J_IN,1)
      CALL LOCATE(X_CAL,N_CAL,X_FI,J_FI,2) 
      N_CMP=J_FI-J_IN+1
C
      WRITE(IUO1,50) N_CMP,X_CAL(J_IN),X_CAL(J_FI)
      WRITE(IUO1,30) STEP_CAL
C
C  Formats
C
  10  FORMAT(//,10X,'<<<<<  N_SIZE SHOULD BE LARGER THAN ',I5,
     1      '  >>>>>',//)
  20  FORMAT(//,6X,A40,'  CONTAINS ',I5,'  POINTS',/,25X,
     1       ' FIRST POINT : ',F8.2,2X,'LAST POINT : ',F8.2)
  30  FORMAT(26X,'STEP :        ',F8.2,//)
  40  FORMAT(26X,'STEP :        IRREGULAR',//)
  50  FORMAT(//,6X,'INTERSECTION FOR THE ANALYSIS         ',
     1       '    CONTAINS ',I5,'  POINTS',/,25X,
     2       ' FIRST POINT : ',F8.2,2X,'LAST POINT : ',F8.2)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE EXPE_INTERPOLATE(X_EXP,I_EXP,X_CAL,I_CAL,N_EXP,N_CMP,
     1                            J_IN,STEP_CAL,IUO1,NUNIT3,
     2                            X2_EXP,I2_EXP)
C
C  This subroutine interpolates the (symmetrized) experimental file
C     on the calculation grid within the intersection bounds
C  The result is stored as (X2_EXP,I2_EXP)
C
C
C  Input parameters:
C
C       X_EXP    :  x coordinates of experimental file
C       I_EXP    :  y coordinates of experimental file 
C       X_CAL    :  x coordinates of calculation file
C       I_CAL    :  y coordinates of calculation file 
C       N_EXP    :  number of points in experimental curve
C       N_CMP    :  number of points in intersection
C       J_IN     :  position of first point of intersection in 
C                   calculation grid
C       STEP_CAL :  step of calculation file
C       IUO1     :  checkfile index for printing
C       NUNIT3   :  index of file where to write the interpolated
C                   experimental points
C
C  Output parameters:
C
C       X2_EXP   :  x coordinates of interpolated experimental file
C       I2_EXP   :  y coordinates of interpolated experimental file 
C                    
C
C
C   Author :  D. Sébilleau
C  
C                                           Last modified :  5 Sep 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 X_CAL(N_SIZE),I_CAL(N_SIZE)
      REAL*4 X_EXP(N_SIZE),I_EXP(N_SIZE)
      REAL*4 X2_EXP(N_SIZE),I2_EXP(N_SIZE)
      REAL*4 YS(N_SIZE)
      REAL*4 YP1,YPN
C
      DATA SMALL /0.0001/
C
C  Initialization
C
      DO J=1,N_SIZE
C
        YS(J)=0.
C
      ENDDO
C
      IF(ABS(I_EXP(1)-I_CAL(J_IN)).GT.SMALL) THEN
C
        DO J=1,N_CMP
          X2_EXP(J)=X_CAL(J_IN)+FLOAT(J-1)*STEP_CAL
        ENDDO
C
        YP1=(-I_EXP(3)+4.*I_EXP(2)-3.*I_EXP(1))/(2.*STEP_CAL)
        YPN=(3.*I_EXP(N_EXP)-4.*I_EXP(N_EXP-1)+
     1       I_EXP(N_EXP-2))/(2.*STEP_CAL)
C
        CALL SPLINE(X_EXP,I_EXP,N_EXP,YP1,YPN,YS)
C
        DO J=1,N_CMP
          XX=X2_EXP(J)
          CALL SPLINT(X_EXP,I_EXP,YS,N_EXP,XX,YY,*5)
          I2_EXP(J)=YY
          WRITE(NUNIT3,*) X2_EXP(J),I2_EXP(J)
        ENDDO
C
      ELSE
C
        DO J=1,N_CMP
          I2_EXP(J)=I_EXP(J)
          X2_EXP(J)=X_EXP(J)
        ENDDO
C
      ENDIF
C
      GOTO 15
C
   5  WRITE(IUO1,10)
      STOP
C
C  Format
C
  10  FORMAT(//,10X,'<<<<<  WRONG VALUE FOR XA IN SUBROUTINE ',
     1       '   SPLINT  >>>>>',//)                                                  
C
  15  RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE SYMMETRIZE(X,Y,N_POINTS,I_SYM,SYM,STEP,IUO1,NUNIT,
     1                      CHFILE)
C
C  This subroutine symmetrizes the curve Y=f(X) about the SYM axis. 
C    The curve is gievn by the arrays (X(I),Y(I)) with I = 1, N_POINTS
C
C  Input parameters:
C
C       * I_SYM   : type of symmetrization
C       * SYM     : symmetry axis (in degrees)
C       * STEP    : x-step of the input curve
C       * IUO1    : output check file number for printing
C       * NUNIT   : unit number to write the symmetrized file
C       * CHFILE  : name of the output symmetrized file
C
C   Author :  D. Sébilleau
C              
C                                          Last modified : 19 Aug 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 X(N_SIZE),Y(N_SIZE)
      REAL*4 SYM,DLOW,DHIGH,XNEW,YNEW
C
      CHARACTER*48 CHFILE
C
      DATA SMALL               /0.0001/
C
C  Symmetrization of experimental file whenever necessary
C     JSYM is the index of the point nearest to the 
C                symmetry axis SYM
C
      IF(I_SYM.EQ.1) THEN
C
        CALL LOCATE(X,N_POINTS,SYM,JSYM,1)
        DLOW=SYM-X(JSYM)
        DHIGH=X(JSYM)-SYM
        IF(DLOW.GT.DHIGH) THEN
          JSYM=JSYM+1
        ENDIF
        N_SYM=0
C
        DO K=JSYM-1,1,-1
            IF((Y(JSYM-K).GT.SMALL).AND.(Y(JSYM+K).GT.SMALL)) THEN
              YNEW=(Y(JSYM-K)+Y(JSYM+K))/2.
            ELSE
              YNEW=Y(JSYM-K)+Y(JSYM+K)
            ENDIF
            Y(JSYM+K)=YNEW
            Y(JSYM-K)=YNEW
            WRITE(NUNIT,*) X(JSYM-K),Y(JSYM-K)
            N_SYM=N_SYM+1
        ENDDO
C
        WRITE(NUNIT,*) X(JSYM),Y(JSYM)
        N_SYM=N_SYM+1
C  
        DO K=1,JSYM-1
          XNEW=X(JSYM)+FLOAT(K)*STEP
          WRITE(NUNIT,*) XNEW,Y(JSYM+K)
          N_SYM=N_SYM+1
        ENDDO
C
        REWIND NUNIT
C
        DO JLINE=1,N_SIZE
          READ(NUNIT,*,END=100) X(JLINE),Y(JLINE)
        ENDDO
        CLOSE(NUNIT)
C
 100    WRITE(IUO1,10) CHFILE,N_SYM,X(1),X(N_SYM)
        N_EXP=N_SYM
C
      ELSEIF(I_SYM.EQ.-1) THEN
C
C   Checking the possibility to do it : need of either 0 
C              or 2 x SYM
C
        ISYM=0
        IF(X(1).LE.0.) THEN
          ISYM=1
        ELSEIF(X(N_EXP).GE.2*SYM) THEN
          ISYM=2
        ENDIF
        IF(ISYM.EQ.0) GOTO 500
C
        CALL LOCATE(X,N_POINTS,SYM,JSYM,1)
        DLOW=SYM-X(JSYM)
        DHIGH=X(JSYM)-SYM
        IF(DLOW.GT.DHIGH) THEN
          JSYM=JSYM+1
        ENDIF
        IF(ISYM.EQ.1) THEN
          CALL LOCATE(X,N_POINTS,0.,JSYM0,2)
          IF(ABS(X(JSYM0)).GT.ABS(X(JSYM0+1))) THEN
            JSYM0=JSYM0+1
          ENDIF
          JDELTA=JSYM-JSYM0
        ELSEIF(ISYM.EQ.2) THEN
          CALL LOCATE(X,N_POINTS,2.*SYM,JSYM1,1)
          JDELTA=JSYM1-JSYM
        ENDIF
C
        N_SYM=0
C
        DO K=1,JSYM
          IF((Y(K).GT.SMALL).AND.(Y(JDELTA+K).GT.SMALL)) THEN
            YNEW=(Y(K)+Y(JDELTA+K))/2.
          ELSE
            YNEW=(Y(K)+Y(JDELTA+K))
          ENDIF
          Y(K)=YNEW
          Y(JDELTA+K)=YNEW
          WRITE(NUNIT,*) X(K),Y(K)
          N_SYM=N_SYM+1
        ENDDO
C
        WRITE(NUNIT,*) X(JSYM),Y(JSYM)
        N_SYM=N_SYM+1
C
        DO K=1,JSYM-2
            J_NEW=JSYM+K
            XNEW=X(JSYM)+FLOAT(K)*STEP
            IF(XNEW.GT.360.) GOTO 200
            WRITE(NUNIT,*) XNEW,Y(J_NEW)
            N_SYM=N_SYM+1
        ENDDO
C
 200    REWIND NUNIT
C
        DO JLINE=1,N_SIZE
          READ(NUNIT,*,END=300) X(JLINE),Y(JLINE)
        ENDDO
        CLOSE(NUNIT)
C
 300    WRITE(IUO1,10) CHFILE,N_SYM,X(1),X(N_SYM)
        N_EXP=N_SYM
C
        GOTO 400
 500    WRITE(IUO1,20)
        I_SYM=0
C
      ENDIF
C
C  Formats
C
  10  FORMAT(6X,'SYMMETRIZED FILE : ',A44,/,48X,'CONTAINS ',I5,
     1       '  POINTS',/,25X,' FIRST POINT : ',F8.2,2X,
     2       'LAST POINT : ',F8.2)
  20  FORMAT(//,10X,'-->  IMPOSSIBLE TO SYMMETRIZE : CARRYING ON ',
     1      'WITHOUT',//)
C
 400  RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE MODULATION_FUNCTION(X,Y,N_CMP,NUNIT,IUO1,NORM,
     1                               I_CAL,I_PRINT)
C
C  This subroutine computes the modulation function corresponding
C    to a given input file
C
C
C  Input parameters:
C
C       X        :  x coordinates of the input file
C       Y        :  y coordinates of the input file
C       N_CMP    :  number of points in intersection
C       IUO1     :  checkfile index for printing
C       NORM     :  switch for normalization
C       I_CAL    :  switch experiment/calculation
C       I_PRINT  :  switch for printing
C
C  Output parameters:
C
C       X        :  x coordinates of the output file
C       Y        :  y coordinates of the output file
C                    
C
C
C   Author :  D. Sébilleau
C  
C                                           Last modified :  9 Sep 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 X(N_SIZE),Y(N_SIZE)
      REAL*4 SUM_Y,AVER_Y
C
      SUM_Y=0.
C
      DO J=1,N_CMP
        SUM_Y=SUM_Y+Y(J)
      ENDDO
C
      AVER_Y=SUM_Y/FLOAT(N_CMP)
C
      IF(I_PRINT.EQ.1) THEN
        IF(I_CAL.EQ.0) THEN
          WRITE(IUO1,10) AVER_Y
        ELSE
          WRITE(IUO1,20) AVER_Y
        ENDIF
      ENDIF
C
      DO J=1,N_CMP
        Y(J)=(Y(J)-ABS(AVER_Y))/ABS(AVER_Y)
        WRITE(NUNIT,*) X(J),Y(J)
      ENDDO
      IF(NORM.GT.0) THEN
        CALL NORMALIZE_CURVE(X,Y,N_CMP,NORM,IUO1)
      ENDIF
C
      CLOSE(NUNIT)
C
C  Formats
C
  10  FORMAT(10X,'----> AVERAGE VALUE FOR EXPERIMENT  :  ',E12.6)
  20  FORMAT(10X,'----> AVERAGE VALUE FOR CALCULATION :  ',E12.6)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE DISTRIBUTION_FUNCTION(X,Y,N_CMP,NUNIT)
C
C  This subroutine computes the cumulative distribution function 
C    corresponding to a given input file
C
C
C  Input parameters:
C
C       X        :  x coordinates of the input file
C       Y        :  y coordinates of the input file
C       N_CMP    :  number of points in file
C
C  Output parameters:
C
C       X        :  x coordinates of the output file
C       Y        :  y coordinates of the output file
C                    
C
C
C   Author :  D. Sébilleau
C  
C                                           Last modified :  5 Sep 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 X(N_SIZE),Y(N_SIZE),F(N_SIZE)
      REAL*4 SUM_Y,CUM_Y
C
C  Initialization
C
      DO J=1,N_CMP
C
        F(J)=0.
C
      ENDDO
C
      SUM_Y=0.
C
      DO J=1,N_CMP
        SUM_Y=SUM_Y+Y(J)
      ENDDO
C
      CUM_Y=0.
C
      DO J=1,N_CMP
        CUM_Y=CUM_Y+Y(J)
        F(J)=CUM_Y/SUM_Y
        WRITE(NUNIT,*) X(J),F(J)
      ENDDO
C
      DO J=1,N_CMP
        Y(J)=F(J)
      ENDDO            
C
      CLOSE(NUNIT)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE CURVATURE_FUNCTION(X,Y,N_POINTS,NUNIT)
C
C  This subroutine computes the curvature function 
C    corresponding to a given input file
C
C
C  Input parameters:
C
C       X        :  x coordinates of the input file
C       Y        :  y coordinates of the input file
C       N_POINTS :  number of points in file
C
C  Output parameters:
C
C       X        :  x coordinates of the output file
C       Y        :  y coordinates of the output file
C                    
C
C
C   Author :  D. Sébilleau
C  
C                                           Last modified : 16 Sep 2015
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 X(N_SIZE),Y(N_SIZE)
      REAL*4 S(N_SIZE),K(N_SIZE)
      REAL*4 T(N_SIZE,2),N(N_SIZE,2)
      REAL*4 S1(N_SIZE),IA(N_SIZE)
C
C  Calling the curve parameters subroutine
C
      CALL CURVE_PARAM(X,Y,N_POINTS,S,S1,K,IA,T,N)
C
      DO J=1,N_POINTS
        Y(J)=K(J)
        WRITE(NUNIT,*) X(J),Y(J)
      ENDDO            
C
      CLOSE(NUNIT)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE CFHWR_FUNCTION(X,Y,N_POINTS,NUNIT)
C
C  This subroutine computes the Cui-Femiani-Hu-Wonka-Razdan signature function 
C    corresponding to a given input file
C
C
C  Input parameters:
C
C       X        :  x coordinates of the input file
C       Y        :  y coordinates of the input file
C       N_POINTS :  number of points in file
C
C  Output parameters:
C
C       X        :  x coordinates of the output file
C       Y        :  y coordinates of the output file
C                    
C
C
C   Author :  D. Sébilleau
C  
C                                           Last modified : 16 Sep 2015
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 X(N_SIZE),Y(N_SIZE)
      REAL*4 S(N_SIZE),K(N_SIZE)
      REAL*4 T(N_SIZE,2),N(N_SIZE,2)
      REAL*4 S1(N_SIZE),IA(N_SIZE)
      REAL*4 IA2(N_SIZE),K2(N_SIZE)
      REAL*4 H
C
C  Calling the curve parameters subroutine
C
      CALL CURVE_PARAM(X,Y,N_POINTS,S,S1,K,IA,T,N)
C
C  Interpolating the Cui-Femiani-Hu-Wonka-Razdan signature function
C
      H=X(2)-X(1)
C
      CALL CFHWR_INTERP(IA,K,S1,N_POINTS,H,IA2,K2)
C
      DO J=1,N_POINTS
        X(J)=IA2(J)
        Y(J)=K2(J)
        WRITE(NUNIT,*) X(J),Y(J)
      ENDDO            
C
      CLOSE(NUNIT)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE CURVE_PARAM(X,Y,N_POINTS,S,S1,K,IA,T,N)
C
C  This subroutine computes the arc length, curvature, integral of  
C    (absolute) curvature, tangent vector and normal vector 
C    of a curve (X,Y)
C
C  It is assumed that the X value are regularly spaced
C
C
C  Input parameters:
C
C       X           :  x coordinates of the input file
C       Y           :  y coordinates of the input file
C       N_POINTS    :  number of points in the file
C
C
C  Output parameters:
C
C       S           :  arc length
C       S1          :  first derivative of arc length
C       K           :  curvature
C       IA          :  integral of |K| over S
C       T           :  tangent vector
C       N           :  normal vector
C
C
C  Reference: M. Cui, J. Femiani, J. Hu, P. Wonka and A. Razdan,
C             Pattern Recognition Letters 30, 1-10 (2009)
C
C   Author :  D. Sébilleau
C
C                                          Last modified : 16 Jan 2015
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 Y(N_SIZE),X(N_SIZE)
      REAL*4 S(N_SIZE),K(N_SIZE)
      REAL*4 T(N_SIZE,2),N(N_SIZE,2)
      REAL*4 TX(N_SIZE),TY(N_SIZE),TX1(N_SIZE),TY1(N_SIZE)
      REAL*4 NX(N_SIZE),NY(N_SIZE),NX1(N_SIZE),NY1(N_SIZE)
      REAL*4 Y1(N_SIZE),Y2(N_SIZE)
      REAL*4 S1(N_SIZE)
      REAL*4 Y3(N_SIZE),Y4(N_SIZE),Y5(N_SIZE)
      REAL*4 F2(N_SIZE),F3(N_SIZE),F4(N_SIZE),F5(N_SIZE)
      REAL*4 F_1(N_SIZE),F_2(N_SIZE),F_3(N_SIZE)
      REAL*4 F_4(N_SIZE),F_5(N_SIZE)
      REAL*4 K_STA(N_SIZE),K_STB(N_SIZE),K_DIR(N_SIZE)
      REAL*4 G(N_SIZE),IA(N_SIZE)
      REAL*4 RES
C
      CHARACTER*3 METH
C
C  Parameters for derivatives calculations
C
      STEP=X(2)-X(1)
      N_CALC=2
      I_FLAG=2
C
C  Parameters for integration calculations
C
      METH='NCQ'
      N_RULE=2
C
C  Initialisations
C
      DO I=1,N_POINTS
C
C.....  Derivatives of input function  .....
C
        Y1(I)=0.0
        Y2(I)=0.0
        Y3(I)=0.0
        Y4(I)=0.0
        Y5(I)=0.0
C
C.....  Arc length and its derivatives and integrals  .....
C
        S(I)=0.0
        S1(I)=0.0
        F_1(I)=0.0
        F_2(I)=0.0
        F_3(I)=0.0
        F_4(I)=0.0
        F_5(I)=0.0
C
C.....  Curvature, tangent vector and normal vector  .....
C
        K(I)=0.0
        TX(I)=0.0
        TY(I)=0.0
        NX(I)=0.0
        NY(I)=0.0
C
        T(I,1)=0.0
        T(I,2)=0.0
        N(I,1)=0.0
        N(I,2)=0.0
C
C.....  Curvature from Serret-Frenet equations  .....
C
        K_DIR(I)=0.0
        K_STA(I)=0.0
        K_STB(I)=0.0
C
C.....  Cui-Femiani-Hu-Wonka-Razdan signature function  .....
C
        G(I)=0.0
        IA(I)=0.0
C
      ENDDO
C
C  Computing Y'(X) and Y''(X) (derivatives with respect to X)
C
      CALL DERIV(Y,N_POINTS,Y1,Y2,Y3,Y4,Y5,N_CALC,STEP,I_FLAG)
C
C  Storage of S1(X), first derivative of arc length S(X)
C
      DO I=1,N_POINTS
C
        S1(I)=SQRT(1.0+Y1(I)*Y1(I))
C
      ENDDO
C
C  Computing the arc length S(X)
C
C
      IF(METH.EQ.'EMS') THEN
C
C.....  Computing the derivatives of S1  .....
C
        I_FL=5
        N_CA=3
        CALL DERIV(S1,N_POINTS,F_1,F_2,F_3,F_4,F_5,N_CA,STEP,I_FL)
C
      ENDIF
C
C.....  Performing the integration  .....
C
      DO I=1,N_POINTS
C
        CALL INTEGR_I(X,S1,F_1,F_3,F_5,1,I,N_POINTS,METH,N_RULE,RES)
        S(I)=RES
C
      ENDDO
C
C
C
C  Computing the curvature K_DIR(X), tangent vector T(X,2)
C             and normal vector N(X,2)
C
      DO I=1,N_POINTS
C
        K_DIR(I)=ABS(Y2(I))/(S1(I)*S1(I)*S1(I))
C
        TX(I)=1.0/S1(I)
        TY(I)=Y1(I)/S1(I)
        NX(I)=-Y1(I)/SIGN(S1(I),Y2(I))
        NY(I)=1.0/SIGN(S1(I),Y2(I))
C
        T(I,1)=TX(I)
        T(I,2)=TY(I)
        N(I,1)=NX(I)
        N(I,2)=NY(I)
C
      ENDDO
C
C  Alternative calculation of curvature through the Serret-Frenet equations:
C
C                  K_STA(X) =  T'.N / S'(X)
C                  K_STB(X) = -N'.T / S'(X)
C
C       F. Mokhtarian and A. K. Mackworth, IEEE Transactions on
C          Pattern Analysis and Machine Intelligence 14, 789 (1992)
C      
      CALL DERIV(TX,N_POINTS,TX1,F2,F3,F4,F5,N_CALC,STEP,I_FLAG)
      CALL DERIV(TY,N_POINTS,TY1,F2,F3,F4,F5,N_CALC,STEP,I_FLAG)
      CALL DERIV(NX,N_POINTS,NX1,F2,F3,F4,F5,N_CALC,STEP,I_FLAG)
      CALL DERIV(NY,N_POINTS,NY1,F2,F3,F4,F5,N_CALC,STEP,I_FLAG)
C
C  We take for K(X) the average of the three ways to compute:
C
      DO I=1,N_POINTS
C
        K_STA(I)=ABS((TX1(I)*NX(I)+TY1(I)*NY(I))/S1(I))
        K_STB(I)=ABS(-(NX1(I)*TX(I)+NY1(I)*TY(I))/S1(I))
C
        K(I)=(K_DIR(I)+K_STA(I)+K_STB(I))/3.0
C
      ENDDO
C
C  Computing the integral of the absolute value of the curvature K 
C            (integration over the arc length S) 
C        Cui-Femiani-Hu-Wonka-Razdan signature function
C
C..... Computing the integrand  .....
C
      DO I=1,N_POINTS
C
        G(I)=ABS(Y2(I))/(S1(I)*S1(I))
C
      ENDDO      
C
C..... Computing the integral
C
      DO I=1,N_POINTS
C
        CALL INTEGR_I(X,G,F_1,F_3,F_5,1,I,N_POINTS,METH,N_RULE,RES)
        IA(I)=RES
C
      ENDDO
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE CFHWR_INTERP(IA,K,S1,N_POINTS,H,IA2,K2)
C
C  This subroutine interpolates the Cui-Femiani-Hu-Wonka-Razdan 
C    signature function (curvature Y = curvature K plotted as a function  
C    of X = int |K| ds, with s the arc length parameter) on a regular grid
C
C  Input parameters:
C
C       IA          :  X as a function of J (= x)
C       K           :  Y as a function of J (= x)
C       S1          :  first derivative S'(x) of arc length parameter S
C       N_POINTS    :  number of points in the file
C       H           :  regular step in x
C
C
C  Output parameters:
C
C       IA2         :  new regularly spaced X
C       K2          :  curvature K on grid X
C
C
C  K2=f(IA2) is the Cui-Femiani-Hu-Wonka-Razdan signature function
C
C
C  Reference: M. Cui, J. Femiani, J. Hu, P. Wonka and A. Razdan,
C             Pattern Recognition Letters 30, 1-10 (2009)
C
C   Author :  D. Sébilleau
C
C                                          Last modified : 16 Jan 2015
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 IA(N_SIZE),IA2(N_SIZE),S1(N_SIZE)
      REAL*4 K(N_SIZE),K2(N_SIZE),KK(N_SIZE)
      REAL*4 X_MIN,X_MAX,STEP
      REAL*4 F1_MIN,F1_MAX,KP_MIN,KP_MAX
C
      DATA IUO1        /6/
C
C  Computing the minimum and maximum of the abscissa IA(J)
C
      X_MIN=1.E+15
      X_MAX=0.0
C
      DO J=1,N_POINTS
C
        X_MIN=MIN(X_MIN,IA(J))
        X_MAX=MAX(X_MAX,IA(J))
C
      ENDDO
C
      STEP=(X_MAX-X_MIN)/FLOAT(N_POINTS-1)
C
C  Dividing the interval into regularly spaced points
C
      DO J=1,N_POINTS
C
        IA2(J)=X_MIN+FLOAT(J-1)*STEP
C
      ENDDO
C
C  Interpolation of curvature K onto this new grid using 
C              cubic spline interpolation
C
C  For this, we need the first derivative dY/dX at the first 
C    point X_MIN at and the last points X_MAX
C
C  We use the fact that dK/dX = dK/dx * dx/dX
C
C  Because int K ds = int K(x) S'(x) dx (K is always >= 0), 
C    we have dX/dx = K(x) S'(x)
C
      F1_MIN=0.5*(-K(3)+4.0*K(2)-3.0*K(1))/H
      F1_MAX=0.5*(3.0*K(N_POINTS)-4.0*K(N_POINTS-1)+K(N_POINTS-2))/H
C
      KP_MIN=F1_MIN/(K(1)*S1(1))
      KP_MAX=F1_MAX/(K(N_POINTS)*S1(N_POINTS))
C
      CALL SPLINE(IA,K,N_POINTS,KP_MIN,KP_MAX,KK)
C
      K2(1)=K(1)
      K2(N_POINTS)=K2(N_POINTS)
C
      DO J=1,N_POINTS
C
        XX=IA2(J)
        CALL SPLINT(IA,K,KK,N_POINTS,XX,YY,*5)
        K2(J)=YY
C
      ENDDO
C
      GOTO 15
C
   5  WRITE(IUO1,10)
      STOP
C
C  Format
C
  10  FORMAT(//,10X,'<<<<<  WRONG VALUE FOR XA IN SUBROUTINE ',
     1       '   SPLINT  >>>>>',//)                                                  
C
C
  15  RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE COMPUTE_SHIFT(EXPE,CALCULATION,N_CMP,NFILE,IUO1,SHIFT)
C
C  This subroutine computes the minimum of experiment 
C     and of all calculations. It will be used to shift 
C     all curves so that all y coordinates are positive 
C     whenever required
C
C
C  Input parameters:
C
C       EXPE        : y coordinates of experimental file
C       CALCULATION : y coordinates of all calculated files
C       N_CMP       : number of points in intersection
C       NFILE       : number of calculation files
C       IUO1        : checkfile index for printing
C
C  Output parameters:
C
C       SHIFT       : shift to be applied to all y coordinates
C                    
C
C
C   Author :  D. Sébilleau
C  
C                                           Last modified :  5 Sep 2014
C
      PARAMETER (N_SIZE=1000,N_FILES=1000)
C
      REAL*4 EXPE(N_SIZE),CALCULATION(N_SIZE,N_FILES)
      REAL*4 MIN_EXP,MIN_CAL,MINIMUM,SHIFT
C
C  Computing the minimum of experiment and of all calculations
C
      SHIFT=0.0
      MIN_EXP=1.0E+30
      DO J=1,N_CMP
        MIN_EXP=MIN(MIN_EXP,EXPE(J))
      ENDDO
C
      MIN_CAL=1.0E+30
      DO JFILE=1,NFILE
        DO J=1,N_CMP
          MIN_CAL=MIN(MIN_CAL,CALCULATION(J,JFILE))
        ENDDO
      ENDDO
C
      MINIMUM=MIN(MIN_EXP,MIN_CAL)
      IF(MINIMUM.LT.0.) THEN
        SHIFT=ABS(MINIMUM)+1.0
        WRITE(IUO1,10) SHIFT
      ENDIF
C
C  Format
C
 10   FORMAT(//,6X,'----> ALL CURVES HAVE BEEN SHIFTED BY ',E12.6,/,15X,
     1       ' TO AVOID NEGATIVE VALUES',//)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE R_FACTOR_1(I_EXP,I_CAL,W,N_POINTS,J_START,CNORM,
     1                      RF1,RF2,RF3,RF4,RF5,I_DEV)
C
C   This subroutine computes R-factors :
C
C         I_DEV = 0    :   5 R-factors without derivatives 
C                                  RF1,RF2,RF3,RF4,RF5
C                          I_EXP and I_CAL are the spectra
C
C         I_DEV = 1    :   2 R-factors with derivatives
C                                  RF3,RF4
C                          I_EXP and I_CAL are the derivatives of the spectra
C
C         I_DEV = 2    :   3 R-factors as defined by D. P. Woodruff et al,
C                          Nucl. Instr. and Meth. in Phys. Res. B 183, 128 (2001)
C
C         SCALEn       :   scaling factor for all R-factors to have the 
C                            same scale (R.J. Koestner, M.A. Van Hove 
C                            and G.A. Somorjai, Surf. Sci. 107, 439 (1981))
C
C  Input parameters:
C
C       I_EXP    :  y coordinates of the experimental file
C       I_CAL    :  y coordinates of the calculation file
C       W        :  weight function
C       N_POINTS :  number of points in the files
C       J_START  :  starting point for calculation (should be 1 !)
C       CNORM    :  scaling coefficient to rescale calculation to experiment
C       I_DEV    :  switch to select the R-factors computed
C
C
C  Output parameters:
C
C       RFn      :  R-factor value
C
C
C   Author :  D. Sébilleau
C
C                                            Last modified : 27 Aug 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 I_EXP(N_SIZE),I_CAL(N_SIZE),W(N_SIZE)
      REAL*4 SUM1,SUM2,SUM3,SUM4,SUMA,SUMB,RF1,RF2,RF3,RF4,CNORM
      REAL*4 SCALE1,SCALE2,SCALE3,SCALE4
C  
C  Scaling factors
C
      SCALE1=0.75
      SCALE2=0.50
      SCALE3=0.75
      SCALE4=0.50  
      SCALE5=1.00 
C
      SUM1=0.  
      SUM2=0. 
      SUM3=0. 
      SUM4=0.
      SUMA=0.
      SUMB=0.
C  
      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
        IF(I_DEV.LE.1) THEN
          SUMA=SUMA+W(J)*ABS(I_EXP(J))
          SUMB=SUMB+W(J)*I_EXP(J)*I_EXP(J)
          SUM1=SUM1+W(J)*ABS(I_EXP(J)-CNORM*I_CAL(JJ))
          SUM2=SUM2+W(J)*ABS(I_EXP(J)-CNORM*I_CAL(JJ))**2
          SUM3=SUM3+W(J)*ABS(I_EXP(J)-CNORM*I_CAL(JJ))**2
          SUM4=SUM4+W(J)*(I_EXP(J)*I_EXP(J)+
     1                    CNORM*CNORM*I_CAL(JJ)*I_CAL(JJ))
        ELSEIF(I_DEV.EQ.2) THEN
          SUM1=SUM1+W(J)*(ABS(I_EXP(J)-CNORM*I_CAL(JJ))**2)/
     1                    ABS(I_EXP(J))
          SUM2=SUM2+W(J)*ABS(I_EXP(J)-CNORM*I_CAL(JJ))**2
          SUM3=SUM3+W(J)*ABS((I_EXP(J)-CNORM*I_CAL(JJ))/I_EXP(J))**2
        ENDIF
C
      ENDDO
C
      IF(I_DEV.EQ.0) THEN
        RF1=SCALE1*SUM1/FLOAT(N_POINTS)
        RF2=SCALE2*SUM2/FLOAT(N_POINTS)
        RF3=SCALE3*SUM1/SUMA
        RF4=SCALE4*SUM2/SUMB
        RF5=SUM3/SUM4
      ELSEIF(I_DEV.EQ.1) THEN
        RF1=0.
        RF2=0.
        RF3=SCALE3*SUM1/SUMA
        RF4=SCALE4*SUM2/SUMB
      ELSEIF(I_DEV.EQ.2) THEN
        RF1=SCALE5*SUM1/FLOAT(N_POINTS)
        RF2=SCALE5*SQRT(SUM2)/FLOAT(N_POINTS)
        RF3=SCALE5*100.*SQRT(SUM3)/FLOAT(N_POINTS)
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE R_FACTOR_2(I_EXP,I_CAL,I_EXP_1,I_CAL_1,I_EXP_2,
     1                      I_CAL_2,W,N_POINTS,J_START,CNORM,V_I,
     2                      RF1,RF2)
C
C   This subroutine computes the Zanazzi-Jona (RF1) 
C               and Pendry (RF2) R-factors
C
C       SCALEn          :  scaling factor for all R-factors to have the 
C                           same scale (R.J. Koestner, M.A. Van Hove 
C                           and G.A. Somorjai, Surf. Sci. 107, 439 (1981)
C
C       MAXE_1          :  maximum absolute value of the experimental
C                           spectrum derivative
C
C  Input parameters:
C
C       I_EXP    :  y coordinates of the experimental file
C       I_CAL    :  y coordinates of the calculation file
C       I_EXP1   :  y coordinates of the first derivative of experimental file
C       I_CAL1   :  y coordinates of the first derivative of the calculation file
C       I_EXP2   :  y coordinates of the second derivative of experimental file
C       I_CAL2   :  y coordinates of the second derivative of the calculation file
C       W        :  weight function
C       N_POINTS :  number of points in the files
C       J_START  :  starting point for calculation (should be 1 !)
C       CNORM    :  scaling coefficient to rescale calculation to experiment
C       V_I      :  imaginary part of the constant potential (eV)
C                     V_I ~ (E_k)**0.3333 is often a good approximation
C
C
C  Output parameters:
C
C       RFn      :  R-factor value
C
C
C   Author :  D. Sébilleau
C
C                                            Last modified : 27 Aug 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 I_EXP(N_SIZE),I_CAL(N_SIZE),I_EXP_1(N_SIZE),I_CAL_1(N_SIZE)
      REAL*4 I_EXP_2(N_SIZE),I_CAL_2(N_SIZE),Y_EXP(N_SIZE),Y_CAL(N_SIZE)
      REAL*4 W(N_SIZE)
      REAL*4 SUM0,SUM1,SUM2,SUM3,CNORM,SCALE1,SCALE2,V_I,MAXE_1
C  
C  Scaling factors
C
      SCALE1=0.50
      SCALE2=0.50
      SCALE3=0.027
C
C 
C  Calculation of the maximum absolute value of 
C          the experimental spectrum derivative
C
C  Calculation of Pendry's Y(J) functions
C
      SUM0=0.
      SUM1=0.
      MAXE_1=0.
C
      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
        SUM1=SUM1+ABS(I_EXP(J))
        MAXE_1=MAX(MAXE_1,ABS(I_EXP_1(J)))
        Y_EXP(J)=I_EXP_1(J)*I_EXP(J)/
     1           (I_EXP(J)*I_EXP(J)+V_I*V_I*I_EXP_1(J)*I_EXP_1(J))
        Y_CAL(JJ)=I_CAL_1(JJ)*I_CAL(JJ)/
     1           (I_CAL(JJ)*I_CAL(JJ)+V_I*V_I*I_CAL_1(JJ)*I_CAL_1(JJ))
C
      ENDDO
C
      SUM0=SCALE3*SUM1
C
      SUM1=0.  
      SUM2=0.
      SUM3=0.
C  
      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
        SUM1=SUM1+W(J)*ABS(I_EXP_2(J)-CNORM*I_CAL_2(JJ))*
     1                 ABS(I_EXP_1(J)-CNORM*I_CAL_1(JJ))/
     2                (ABS(I_EXP_1(J)+MAXE_1))
        SUM2=SUM2+W(J)*((Y_EXP(J)-Y_CAL(JJ))**2)
        SUM3=SUM3+W(J)*(Y_EXP(J)**2+Y_CAL(JJ)**2)
C
      ENDDO
C
       RF1=SCALE1*SUM1/SUM0  
       RF2=SCALE2*SUM2/SUM3
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE SIM_INDEX(I_EXP,I_CAL,W,N_POINTS,J_START,CNORM,
     1                     SI1,SI2,SI3,SI4,SI5,SI6,SI7,SI8,SI9,SI10,
     2                     SI11,SI12,ALPHA,BETA,N_BIN)
C
C   This subroutine computes similarity indices
C
C     VAR_FUN is the variance of function FUN
C     COVAR is the covariance
C     CNORM is a renormalization coefficient so that CNORM*I_CAL 
C       has the same order of magnitude as I_EXP
C
C     LUM : luminance
C     CON : contrast
C     STR : structure
C        (see R. Dosselmann and D. Y. Xue, Signal, Image and Video
C             Processing, 5, 81 (2011))
C
C   Note : SI3 = 1/n * sum_i [ noise(i) ]
C
C     N_BIN : number of points in a bin
C     NBINS : number of bins
C
C    SI1  : Spectral contrast angle (or cosine similarity, or Carbo index)
C    SI2  : Linear correlation
C    SI3  : Integrated similarity index
C    SI4  : Structural similarity index
C    SI5  : Tversky/Tanimoto index
C    SI6  : Dice index
C    SI7  : Sorensen index
C    SI8  : Czekanowski index
C    SI9  : Weighted cross correlation
C    SI10 : Weighted cross correlation
C    SI11 : Horn coefficient
C    SI12 : Binning similarity index
C
C  Input parameters:
C
C       I_EXP    :  y coordinates of the experimental file
C       I_CAL    :  y coordinates of the calculation file
C       W        :  weight function
C       N_POINTS :  number of points in the files
C       J_START  :  starting point for calculation (should be 1 !)
C       CNORM    :  scaling coefficient to rescale calculation to experiment
C       ALPHA    :  alpha parameter in Tversky/Tanimoto index
C       BETA     :  beta parameter in Tversky/Tanimoto index 
C       N_BIN    :  number of points in a bin
C
C
C  Output parameters:
C
C       SIn      :  similarity index value
C
C
C   Author :  D. Sébilleau
C
C                                          Last modified : 23 Sep 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 I_EXP(N_SIZE),I_CAL(N_SIZE),I2_EXP(N_SIZE),I2_CAL(N_SIZE)
      REAL*4 W(N_SIZE)
      REAL*4 VAR_EXP,VAR_CAL,COVAR
      REAL*4 LUM,CON,STR,C1,C2,C3
      REAL*4 X1,X2,X3
      REAL*4 SUM1,SUM2,SUM3,SUM4,SUM5
      REAL*4 SUM_EC_1,SUM_EE_1,SUM_CC_1
      REAL*4 SUM_EC_2,SUM_EE_2,SUM_CC_2
      REAL*4 SUM9_1,SUM10_1,SUM9_2,SUM10_2,SUM9_3,SUM10_3
      REAL*4 SUM_EXP,SUM_CAL,SCAL
      REAL*4 MEAN_EXP,MEAN_CAL,CNORM
      REAL*4 H_EXP,H_CAL,TMP_EXP,TMP_CAL,H_MAX,H_MIN,H_OBS,SUM_EC,TMP_EC
      REAL*4 SUM_BIN_EXP,SUM_BIN_CAL,S_EXPCAL,SUM_SI,SI_N
      REAL*4 SI1,SI2,SI3,SI4,SI5,SI6,SI7,SI8,SI9,SI10,SI11,SI12
C
      INTEGER I_START(N_SIZE),I_END(N_SIZE)
C
      SCAL=1./FLOAT(N_POINTS)
      C1=2.0
      C2=3.0
      C3=0.5*C2
C
C  Defining partial sums, full sums, mean and overall minimum
C
      SUM1=0.
      SUM2=0.
      SUM_EXP=0.
      SUM_CAL=0.
C
      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
        SUM_EXP=SUM_EXP+W(J)*I_EXP(J)
        SUM_CAL=SUM_CAL+W(J)*I_CAL(JJ)
        I2_EXP(J)=SUM_EXP
        I2_CAL(JJ)=SUM_CAL
C
      ENDDO
C
      MEAN_EXP=SCAL*SUM_EXP
      MEAN_CAL=SCAL*SUM_CAL
C
C  Computing H(S) Shannon entropy for Horn index (independent of CNORM)
C
      H_EXP=0.
      H_CAL=0.
      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
        TMP_EXP=SUM_EXP/I_EXP(J)
        TMP_CAL=SUM_CAL/I_CAL(JJ)
        H_EXP=H_EXP+W(J)*LOG(TMP_EXP)/TMP_EXP
        H_CAL=H_CAL+W(J)*LOG(TMP_CAL)/TMP_CAL
C
      ENDDO
C
C  Computing H_max, H_min and H_obs
C
      SUM_EC=SUM_EXP+CNORM*SUM_CAL
      H_MAX=0.
      H_OBS=0.
      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
        TMP_EC=I_EXP(J)+CNORM*I_CAL(JJ)
        H_MAX=H_MAX+W(J)*I_EXP(J)*LOG(SUM_EC/I_EXP(J))/SUM_EC+
     1              W(J)*CNORM*I_CAL(JJ)*LOG(SUM_EC/I_CAL(JJ))/SUM_EC
        H_OBS=H_OBS+TMP_EC*LOG(SUM_EC/TMP_EC)/SUM_EC
C
      ENDDO
      H_MIN=SUM_EXP*H_EXP/SUM_EC + CNORM*SUM_CAL*H_CAL/SUM_EC
      SI11=(H_MAX-H_OBS)/(H_MAX-H_MIN)
C
C  Binning method:
C
C  Initial and final value of each of the NBINS bins
C
      NBINS=N_POINTS/N_BIN
C
      I_START(1)=1
      I_END(1)=N_BIN
C
      SUM_SI=0.
C
      DO N=2,NBINS
C
        I_START(N)=I_END(N-1)+1
        IF(N.LT.NBINS) THEN
          I_END(N)=I_START(N)+N_BIN-1
        ELSE
          I_END(N)=N_POINTS
        ENDIF
C
      ENDDO
C
C  Computation of S_EXPCAL
C
      DO N=1,NBINS
C
        S_EXPCAL=0.
C
        DO J=1,N
C
          SUM_BIN_EXP=0.
          SUM_BIN_CAL=0.
C
C    Sums within bin J
C
          DO I=I_START(J),I_END(J)
C
            SUM_BIN_EXP=SUM_BIN_EXP+I_EXP(I)
            SUM_BIN_CAL=SUM_BIN_CAL+CNORM*I_CAL(I)
C
          ENDDO
C
          S_EXPCAL=S_EXPCAL+MIN(SUM_BIN_EXP,SUM_BIN_CAL)
C
        ENDDO
C
C
        SI_N=S_EXPCAL/(SUM_EXP+CNORM*SUM_CAL-S_EXPCAL)
        SUM_SI=SUM_SI+SI_N
C
      ENDDO
C
      SI12=SUM_SI/FLOAT(NBINS)
C
      SUM_EC_1=0.
      SUM_EE_1=0.
      SUM_CC_1=0.
      SUM_EC_2=0.
      SUM_EE_2=0.
      SUM_CC_2=0.
C
C  Cross correlation coefficients
C
C
      DO I=1,N_POINTS
C
        SUM9_1=0.
        SUM10_1=0.
        SUM9_2=0.
        SUM10_2=0.
        SUM9_3=0.
        SUM10_3=0.
C
        DO J=1,N_POINTS
C
          JJ=J_START+J-1
C
          IF(J.LE.N_POINTS-I) THEN
            SUM9_1=SUM9_1+I_EXP(J)*I_CAL(JJ+I)
            SUM9_2=SUM9_2+I_EXP(J)*I_EXP(J+I)
            SUM9_3=SUM9_3+I_CAL(JJ)*I_CAL(JJ+I)
            SUM10_1=SUM10_1+I_EXP(J)*I_CAL(JJ+I)
            SUM10_2=SUM10_2+I_EXP(J)*I_EXP(J+I)
            SUM10_3=SUM10_3+I_CAL(JJ)*I_CAL(JJ+I)
          ELSE
            SUM10_1=SUM10_1+I_EXP(J)*I_CAL(JJ+I-N_POINTS)
            SUM10_2=SUM10_2+I_EXP(J)*I_EXP(J+I-N_POINTS)
            SUM10_3=SUM10_3+I_CAL(JJ)*I_CAL(JJ+I-N_POINTS)
          ENDIF
C
        ENDDO
C
        SUM_EC_1=SUM_EC_1+W(I)*SUM9_1
        SUM_EE_1=SUM_EE_1+W(I)*SUM9_2
        SUM_CC_1=SUM_CC_1+W(I)*SUM9_3
        SUM_EC_2=SUM_EC_2+W(I)*SUM10_1
        SUM_EE_2=SUM_EE_2+W(I)*SUM10_2
        SUM_CC_2=SUM_CC_2+W(I)*SUM10_3
C
      ENDDO
C
      SI9=SUM_EC_1/SQRT(SUM_EE_1*SUM_CC_1)
      SI10=SUM_EC_2/SQRT(SUM_EE_2*SUM_CC_2)
C
C  Defining variance and covariance
C
      VAR_EXP=0.
      VAR_CAL=0.
      COVAR=0.
C
      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
C  Variance and covariance
C
        X1=W(J)*(I_EXP(J)-MEAN_EXP)*(I_EXP(J)-MEAN_EXP)
        X2=W(J)*(I_CAL(JJ)-MEAN_CAL)*(I_CAL(JJ)-MEAN_CAL)
        X3=W(J)*(I_EXP(J)-MEAN_EXP)*(I_CAL(JJ)-MEAN_CAL)
        VAR_EXP=VAR_EXP+X1
        VAR_CAL=VAR_CAL+X2
        COVAR=COVAR+X3
C
      ENDDO
C
      VAR_EXP=SCAL*VAR_EXP
      VAR_CAL=SCAL*CNORM*CNORM*VAR_CAL
      COVAR=SCAL*CNORM*COVAR
C
      LUM=(2.*MEAN_EXP*CNORM*MEAN_CAL+C1)/(MEAN_EXP*MEAN_EXP+
     1     CNORM*CNORM*MEAN_CAL*MEAN_CAL+C1)
      CON=(2.*SQRT(VAR_EXP*VAR_CAL)+C2)/(VAR_EXP+VAR_CAL+C2)
      STR=(COVAR+C3)/(SQRT(VAR_EXP*VAR_CAL)+C3)
C
      SUM1=0.
      SUM2=0.
      SUM3=0.
      SUM4=0.
      SUM5=0.
C
      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
        SUM1=SUM1+W(J)*I_EXP(J)*CNORM*I_CAL(JJ)
        SUM2=SUM2+W(J)*I_EXP(J)*I_EXP(J)
        SUM3=SUM3+W(J)*CNORM*CNORM*I_CAL(JJ)*I_CAL(JJ)
        SUM4=SUM4+W(J)*MIN(I_EXP(J),CNORM*I_CAL(JJ))
        SUM5=SUM5+W(J)*ABS((I2_EXP(J)/SUM_EXP)-(I2_CAL(JJ)/SUM_CAL))
C
      ENDDO
C
      SI1=SCAL*SUM1/SQRT(SUM2*SUM3)
      SI2=SCAL*COVAR/SQRT(VAR_EXP*VAR_CAL)
      SI3=1.0-SCAL*SUM5
      SI4=SCAL*SCAL*LUM*CON*STR
      SI5=SUM1/(ALPHA*SUM2+BETA*SUM3+(1.0-ALPHA-BETA)*SUM1)
      SI6=2.0*SUM1/(SUM2+SUM3)
      SI7=2.0*SUM1/(SUM2+SUM3+2.0*SUM1)
      SI8=2.0*SUM4/(SUM2+SUM3)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE DISTANCE(X_EXP,I_EXP,I_CAL,W,N_POINTS,J_START,CNORM,
     1                    DI1,DI2,DI3,DI4,DI5,DI6,DI7,DI8,DI9,DI10,DI11,
     2                    DI12,DI13,DI14,DI15,DI16,DI17,DI18,DI19,DI20,
     3                    DI21,DI22,DI23,DI24,ALPHA,I_BETA,SIGMA,L)
C
C   This subroutine computes distances between two data sets
C
C     DI1  :  Euclidian
C     DI2  :  Minkowski
C     DI3  :  Taxicab
C     DI4  :  Average linkage
C     DI5  :  Hellinger
C     DI6  :  Kullback-Leibler
C     DI7  :  Mahalanobis
C     DI8  :  Folded spectra
C     DI9  :  Triangle weight
C     DI10 :  Quadratic weight
C     DI11 :  Gaussian weight
C     DI12 :  Hausdorff
C     DI13 :  Levy
C     DI14 :  Battacharyya
C     DI15 :  Canberra
C     DI16 :  Jeffrey
C     DI17 :  Jensen
C     DI18 :  Histogram intersection
C     DI19 :  Soergel
C     DI20 :  Taneja
C     DI21 :  Kumar-Johnson
C     DI22 :  Jensen difference
C     DI23 :  Min-symmetric chi2
C     DI24 :  discrete Frechet distance
C
C  Input parameters:
C
C       X_EXP    :  x coordinates of the experimental file
C       I_EXP    :  y coordinates of the experimental file
C       I_CAL    :  y coordinates of the calculation file
C       W        :  weight function
C       N_POINTS :  number of points in the files
C       J_START  :  starting point for calculation (should be 1 !)
C       CNORM    :  scaling coefficient to rescale calculation to experiment
C       ALPHA    :  alpha parameter in folded spectra distance (DI8)
C       I_BETA   :  beta parameter in folded spectra distance (DI8)
C       SIGMA    :  sigma parameter in Gaussian weight distance (DI11)
C       L        :  l parameter in triangle weight distance (DI9)
C
C
C  Output parameters:
C
C       DIn      :  distance value
C
C
C   Author :  D. Sébilleau
C
C                                          Last modified :  1 Sep 2014
C
      PARAMETER (N_SIZE=1000,N_TRIES=10000)
C
      REAL*4 X_EXP(N_SIZE),I_EXP(N_SIZE),I_CAL(N_SIZE)
      REAL*4 W(N_SIZE),DF(N_SIZE,N_SIZE)
      REAL*4 SUM1,SUM2,SUM3,SUM4,SUM4_1,SUM5,SUM6,SUM7,SUM8,SUM9
      REAL*4 DI1,DI2,DI3,DI4,DI5,DI6,DI7,DI8,DI9,DI10,DI11,DI12
      REAL*4 DI13,DI14,DI15,DI16,DI17,DI18,DI19,DI20,DI21,DI22,DI23,DI24
      REAL*4 DXY,DYX,PXY,WIJ_1,WIJ_2,WIJ_3,WIJ_4,M_IJ,SUM10,SUM11,MI
      REAL*4 SUM14,SUM15,SUM16,SUM17_1,SUM17_2,SUM17_3,SUM20,SUM21,SUM22
      REAL*4 SUM23_1,SUM23_2
      REAL*4 MINJ,MAXJ,MINI,MAXI,MIN123,DIJ,TVD,XMH,XPH,LOW,HIGH,CAL
      REAL*4 SUM_EXP,SUM_CAL,MEAN_EXP,MEAN_CAL,SIGMA
      REAL*4 NORM_EXP,NORM_CAL,MIN_EC,MAX_EC,X_POINTS
C
      COMPLEX*16 COV_MAT(N_SIZE,N_SIZE),WORK(4*N_SIZE)
C
      INTEGER IPIV(N_SIZE)
C
      DATA SMALL,PI /0.001,3.141593/ 
C
      I_MINK=4
      SUM_EXP=0.
      SUM_CAL=0.
      NORM_EXP=0.
      NORM_CAL=0.
      TVD=0.
      X_POINTS=FLOAT(N_POINTS)
C
C  Mean
C
      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
        SUM_EXP=SUM_EXP+I_EXP(J)
        SUM_CAL=SUM_CAL+I_CAL(JJ)
        NORM_EXP=NORM_EXP+I_EXP(J)*I_EXP(J)
        NORM_CAL=NORM_CAL+I_CAL(JJ)*I_CAL(JJ)
C
      ENDDO
C
      NORM_EXP=SQRT(NORM_EXP)
      NORM_CAL=SQRT(NORM_CAL)
      MEAN_EXP=SUM_EXP/X_POINTS
      MEAN_CAL=SUM_CAL/X_POINTS
C
      SUM1=0.
      SUM2=0.
      SUM3=0.
      SUM5=0.
      SUM6=0.
      SUM14=0.
      SUM15=0.
      SUM16=0.
      SUM17_1=0.
      SUM17_2=0.
      SUM17_3=0.
      SUM20=0.
      SUM21=0.
      SUM22=0.
      SUM23_1=0.
      SUM23_2=0.
C
      TVD=0.
      M_IJ=0.
      MAX_EC=0.
      MIN_EC=0.
C
C  Pointwise distances plus total variation distance (TVD)
C    used as un upper bound for Levy distance
C
      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
        DXY=I_EXP(J)-CNORM*I_CAL(JJ)
        PXY=I_EXP(J)*CNORM*I_CAL(JJ)
        MI=0.5*(I_EXP(J)+CNORM*I_CAL(JJ))
        TVD=MAX(TVD,ABS(DXY))
C
        SUM1=SUM1+W(J)*(DXY*DXY)
        SUM2=SUM2+W(J)*(DXY**I_MINK)
        SUM3=SUM3+W(J)*ABS(DXY)
        SUM5=SUM5+W(J)*(SQRT(I_EXP(J))-SQRT(CNORM*I_CAL(JJ)))**2
        SUM6=SUM6+W(J)*I_EXP(J)*LOG(I_EXP(J)/(CNORM*I_CAL(JJ)))
        SUM14=SUM14+W(J)*SQRT(I_EXP(J)*CNORM*I_CAL(JJ))
        SUM15=SUM15+W(J)*ABS(DXY/(MI+MI)) 
        SUM16=SUM16+W(J)*I_EXP(J)*LOG(I_EXP(J)/MI)+CNORM*I_CAL(JJ)*
     1              LOG(CNORM*I_CAL(JJ)/MI)
        SUM17_1=SUM17_1+MI**ALPHA
        SUM17_2=SUM17_2+I_EXP(J)**ALPHA
        SUM17_3=SUM17_3+SQRT(MI*I_EXP(J))**ALPHA
        SUM20=SUM20+W(J)*MI*LOG(MI/SQRT(PXY))
        SUM21=SUM21+W(J)*2.0*DXY*MI*DXY*MI/(PXY**1.5)
        SUM22=SUM22+W(J)*(0.5*(I_EXP(J)*LOG(I_EXP(J))+
     1                         CNORM*I_CAL(JJ)*LOG(CNORM*I_CAL(JJ)))-
     2                         MI*LOG(MI))
        SUM23_1=SUM23_1+DXY*DXY/I_EXP(J)
        SUM23_2=SUM23_2+DXY*DXY/(CNORM*I_CAL(JJ))
C
        MIN_EC=MIN_EC+W(J)*MIN(I_EXP(J),CNORM*I_CAL(JJ))
        MAX_EC=MAX_EC+W(J)*MAX(I_EXP(J),CNORM*I_CAL(JJ))
C
        DO I=1,N_POINTS
C
          II=J_START+I-1
C
          M_IJ=MAX(M_IJ,ABS(I_EXP(J)-CNORM*I_CAL(II)))
C
        ENDDO
C
      ENDDO
C
      DI1=SQRT(SUM1)
      DI2=SUM2**(1./FLOAT(I_MINK))
      DI3=SUM3
      DI5=SQRT(2.0*SUM5)
      DI6=SUM6
      DI14=-ALOG(SUM14)
      DI15=SUM15
      DI16=SUM16
      DI17=LOG(SUM17_3)-0.5*(LOG(SUM17_1)+LOG(SUM17_2))
      DI17=DI17/(1.-ALPHA)
      DI18=1.-MIN_EC/MIN(NORM_EXP,CNORM*NORM_CAL)
      DI19=SUM3/MAX_EC
      DI20=SUM20
      DI21=SUM21
      DI22=SUM22
      DI23=MIN(SUM23_1,SUM23_2)
C
C  Min and Max for Hausdorff distance
C
      MAXI=0.0
      DO I=1,N_POINTS
C
        MINJ=9.E+30
        DO J=1,N_POINTS
C
          JJ=J_START+J-1
C
          DIJ=ABS(I_EXP(I)-CNORM*I_CAL(JJ))
          MINJ=MIN(MINJ,DIJ)
C
        ENDDO
C
        MAXI=MAX(MAXI,MINJ)
C
      ENDDO
C
      MAXJ=0.0
      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
        MINI=9.E+30
        DO I=1,N_POINTS
C
          DIJ=ABS(I_EXP(I)-CNORM*I_CAL(JJ))
          MINI=MIN(MINI,DIJ)
C
        ENDDO
C
        MAXJ=MAX(MAXJ,MINI)
C
      ENDDO
C
      DI12=MAX(MAXI,MAXJ)
C
      SUM4=0.
      SUM8=0.
      SUM9=0.
      SUM10=0.
      SUM11=0.
C
C  Neighbourhood distances and storage of covariance matrix
C
      DO I=1,N_POINTS
C
        II=J_START+I-1
C
        DXY=I_EXP(I)-CNORM*I_CAL(II)
C
        SUM4_1=0.0
        DO J=1,N_POINTS
C
          JJ=J_START+J-1
C
          IMJ=ABS(I-J)
          DYX=I_EXP(J)-CNORM*I_CAL(JJ)
          DIJ=ABS(I_EXP(I)-CNORM*I_CAL(JJ))
          WIJ_1=1.-DIJ/M_IJ
          WIJ_2=1./(1.+ALPHA*FLOAT(IMJ**I_BETA))
          IF(IMJ.LT.L) THEN
            WIJ_3=1.-FLOAT(IMJ)/FLOAT(L)
          ELSE
            WIJ_3=0.
          ENDIF
          WIJ_4=EXP(-DIJ*DIJ/(2.*SIGMA*SIGMA))/(2*PI*SIGMA*SIGMA)
C
          SUM4_1=SUM4_1+ABS(I_EXP(I)-CNORM*I_CAL(JJ))
          SUM8=SUM8+DXY*WIJ_2*DYX
          SUM9=SUM9+DXY*WIJ_3*DYX
          SUM10=SUM10+DXY*WIJ_1*DYX
          SUM11=SUM11+DXY*WIJ_4*DYX
C
          COV_MAT(I,J)=DCMPLX((I_EXP(I)-MEAN_EXP)*
     1                        CNORM*(I_CAL(JJ)-MEAN_CAL))
C
        ENDDO
C
        SUM4=SUM4+W(I)*SUM4_1
C
      ENDDO
C
      DI4=SUM4/(X_POINTS*X_POINTS)
      DI8=SUM8
      DI9=SUM9
      DI10=SUM10
      DI11=SUM11
C
C  Computing Mahalanobis distance using Lapack inversion
C
      LWORK=N_POINTS
C
      CALL ZGETRF(N_POINTS,N_POINTS,COV_MAT,N_SIZE,IPIV,INFO1)
      IF(INFO1.NE.0) THEN
         WRITE(6,*) '     --->  INFO1 =',INFO1
      ELSE
        CALL ZGETRI(N_POINTS,COV_MAT,N_SIZE,IPIV,WORK,LWORK,INFO)
        IF(INFO.NE.0) THEN
          WRITE(6,*) '    --->  WORK(1),INFO =',WORK(1),INFO
        ENDIF
      ENDIF
C
      SUM7=0.
C
      DO I=1,N_POINTS
C
        II=J_START+I-1
C
        DXY=I_EXP(I)-CNORM*I_CAL(II)
C
        DO J=1,N_POINTS
C
          JJ=J_START+J-1
C
           DYX=I_EXP(J)-CNORM*I_CAL(JJ)
C
           SUM7=SUM7+DXY*REAL(REAL(COV_MAT(I,J)))*DYX
C
        ENDDO
C
      ENDDO
C
      DI7=SUM7
C
C  Computing Levy distance starting from TVD and 
C    using a divide and conquer algorithm
C
      DIST0=TVD
      I_GOOD=1
C
      DO J_TRY=1,N_TRIES
C
        IF(I_GOOD.EQ.1) THEN
          DIST=DIST0/2.
        ELSE
          DIST=3.*DIST0/2.
        ENDIF
C
        I_GOOD=1
C
        DO J=1,N_POINTS
C
          JJ=J_START+J-1
C
C  Calculation of I_EXP at x = XMH and XPH :
C    We use a 3-point Lagrange interpolation
C               with step h/2
C
          XMH=X_EXP(I)-DIST
          XPH=X_EXP(I)+DIST
C
          CALL LOCATE(X_EXP,N_POINTS,XMH,JMH,1)
          CALL LOCATE(X_EXP,N_POINTS,XPH,JPH,1)
C
          IF(JMH.EQ.1) JMH=2
          IF(JPH.EQ.1) JPH=2
          IF(JMH.EQ.N_POINTS) JMH=N_POINTS-1
          IF(JPH.EQ.N_POINTS) JPH=N_POINTS-1
C
          RES_MH=(-I_EXP(JMH-1)+6.*I_EXP(JMH)+3.*I_EXP(JMH+1))/8.
          RES_PH=(-I_EXP(JPH-1)+6.*I_EXP(JPH)+3.*I_EXP(JPH+1))/8.
          LOW=RES_MH-DIST
          HIGH=RES_PH+DIST
          CAL=CNORM*I_CAL(JJ)
C
          IF((LOW.GT.CAL).OR.(CAL.GT.HIGH)) THEN
            I_GOOD=0
            GOTO 10
          ENDIF
          
C
        ENDDO
C
        IF((I_GOOD.EQ.1).AND.(ABS(DIST-DIST0).LE.SMALL)) GOTO  20
C
  10    DIST0=DIST
C
      ENDDO
C
  20  DI13=DIST
C
C  Computing discrete Frechet distance using Mosig-Clausen algorithm
C
      DF(1,1)=ABS(I_EXP(1)-CNORM*I_CAL(1))
      DO J=2,N_POINTS
C
        DF(1,J)=MAX(ABS(I_EXP(1)-CNORM*I_CAL(J)),DF(1,J-1))
C
      ENDDO
C
      DO I=2,N_POINTS
C
        DF(I,1)=MAX(ABS(I_EXP(I)-CNORM*I_CAL(1)),DF(I-1,1))
C
        DO J=2,N_POINTS
C
          MIN123=MIN(DF(I,J-1),DF(I-1,J),DF(I-1,J-1))
          DF(I,J)=MAX(ABS(I_EXP(I)-CNORM*I_CAL(J)),MIN123)
C
        ENDDO
C
      ENDDO
C
      DI24=DF(N_POINTS,N_POINTS)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE GOODNESS(I_EXP,I_CAL,W,N_POINTS,J_START,CNORM,
     1                    GF1,GF2,GF3,GF4,GF5,GF6,GF7,GF8,GF9,GF10,GF11,
     1                    GF12,N_BIN,ALPHA)
C
C   This subroutine computes goodness of fit between two data sets
C
C     GF1  :  Pearson chi^2
C     GF2  :  Kolmogorov-Smirnov
C     GF3  :  Kuiper
C     GF4  :  Cramér-Von Mises
C     GF5  :  Anderson-Darling
C     GF6  :  Watson
C     GF7  :  Likelihood ratio
C     GF8  :  Power divergence
C     GF9  :  Freeman-Tukey
C     GF10 :  Cowell
C     GF11 :  Cressie-Read divergence
C     GF12 :  Phi-divergence
C
C  Input parameters:
C
C       I_EXP    :  y coordinates of the experimental file
C       I_CAL    :  y coordinates of the calculation file
C       W        :  weight function
C       N_POINTS :  number of points in the files
C       J_START  :  starting point for calculation (should be 1 !)
C       CNORM    :  scaling coefficient to rescale calculation to experiment
C       N_BIN    :  number of bins in Freeman-Tukey formula (GF9)
C       ALPHA    :  alpha parameter in power divergence (GF8), Cowell (GF10)
C                      Cressie-Read divergence (GF11) and phi-divergence (GF12)
C
C
C  Output parameters:
C
C       GFn      :  goodness of fit value
C
C
C   Author :  D. Sébilleau
C
C                                          Last modified :  1 Sep 2014
C
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 I_EXP(N_SIZE),I_CAL(N_SIZE)
      REAL*4 F_EXP(N_SIZE),F_CAL(N_SIZE),W(N_SIZE)
      REAL*4 CNORM,SUM_EXP,SUM_CAL,CUM_EXP,CUM_CAL,SUM_FDIF
      REAL*4 SUM1,MAX_1,MIN_1,SUM4,SUM5,SUM6,SUM7,SUM8,SUM9,SUM10
      REAL*4 SUM11,SUM12,SUMB,RATIO,RATIO2,MI,F_DIF
      REAL*4 X_POINTS,SUMF_EXP,SUMF_CAL,AVE_EXP,AVE_CAL,EPS
      REAL*4 PHI1,PHI2
C
      DATA EPS /0.0001/
C
      X_POINTS=FLOAT(N_POINTS)
C
C  Sum of experimental/calculation points
C
      CUM_EXP=0.
      CUM_CAL=0.
C

      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
        CUM_EXP=CUM_EXP+I_EXP(J)
        CUM_CAL=CUM_CAL+I_CAL(JJ)
C
      ENDDO
C
C    Calculation of empirical cumulative distribution functions
C
      SUM_EXP=0.
      SUM_CAL=0.
      SUM_FDIF=0.
C
      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
        SUM_EXP=SUM_EXP+I_EXP(J)
        SUM_CAL=SUM_CAL+I_CAL(JJ)
C
        F_EXP(J)=SUM_EXP/CUM_EXP
        F_CAL(J)=SUM_CAL/CUM_CAL
        SUM_FDIF=SUM_FDIF+F_EXP(J)-F_CAL(J)
C
      ENDDO
C
C    Average empirical cumulative distribution functions
C
      SUMF_EXP=0.
      SUMF_CAL=0.
      DO J=1,N_POINTS
C
        SUMF_EXP=SUMF_EXP+F_EXP(J)     
        SUMF_CAL=SUMF_CAL+F_CAL(J)     
C
      ENDDO
C
      AVE_EXP=SUMF_EXP/X_POINTS
      AVE_CAL=SUMF_CAL/X_POINTS
C
C   Calculation of goodness of fit
C
      SUM1=0.
      MAX_1=0.
      MIN_1=9.E+30
      SUM4=0.
      SUM5=0.
      SUM6=0.
      SUM7=0.
      SUM8=0.
      SUM9=0.
      SUM10=0.
      SUM11=0.
      SUM12=0.
C
      DO J=1,N_POINTS
C
        JJ=J_START+J-1
C
        MI=0.5*(I_EXP(J)-CNORM*I_CAL(JJ))
C
        SUM1=SUM1+W(J)*(I_EXP(J)-MI)*(I_EXP(J)-MI)/MI
        F_DIF=F_EXP(J)-F_CAL(J)
        RATIO=F_EXP(J)/F_CAL(J)
        IF(J.LT.N_POINTS) THEN
          RATIO2=(1.-F_EXP(J))/(1.-F_CAL(J))
        ELSE
          RATIO2=1.
        ENDIF
C
        MAX_1=MAX(MAX_1,ABS(F_DIF))
        MIN_1=MIN(MIN_1,ABS(F_DIF))
C
        SUM4=SUM4+W(J)*F_DIF*F_DIF
        IF(J.LT.N_POINTS) THEN
          SUM5=SUM5+W(J)*F_DIF*F_DIF/(F_EXP(J)*(1.-F_EXP(J)))
        ELSE
          SUM5=SUM5+W(J)
        ENDIF
        SUM6=SUM6+W(J)*(F_DIF-SUM_FDIF)*(F_DIF-SUM_FDIF)
        SUM7=SUM7+W(J)*F_EXP(J)*LOG(RATIO)
        SUM8=SUM8+W(J)*F_EXP(J)*((RATIO**ALPHA)-1.)
C
        SUMB=0.
        DO I=1,N_BIN
          SUMB=SUMB+RATIO**(0.5*FLOAT(I))
        ENDDO
C
        SUM9=SUM9+W(J)*SUMB*(SQRT(F_EXP(J))-SQRT(F_CAL(J)))
        SUM10=SUM10+W(J)*(((F_EXP(J)/AVE_EXP)**ALPHA)*
     1                    ((F_CAL(J)/AVE_CAL)**(1.-ALPHA))-1.0)
        SUM11=SUM11+W(J)*(F_EXP(J)*(RATIO**ALPHA)+(1.0-F_EXP(J))*
     1                    (RATIO2**ALPHA)-1.0)
        IF(ABS(ALPHA).LT.EPS) THEN
          PHI1=RATIO-LOG(RATIO)-1.0
          PHI2=RATIO2-LOG(RATIO2)-1.0
        ELSEIF(ABS(ALPHA-1.0).LT.EPS) THEN
          PHI1=RATIO*LOG(RATIO)-RATIO+1.0
          PHI2=RATIO2*LOG(RATIO2)-RATIO2+1.0
        ELSE
          PHI1=(1.-ALPHA+ALPHA*RATIO-(RATIO**ALPHA))/(ALPHA*(1.-ALPHA))
          PHI2=(1.-ALPHA+ALPHA*RATIO2-(RATIO2**ALPHA))/
     1          (ALPHA*(1.-ALPHA))
        ENDIF
        SUM12=SUM12+W(J)*(F_CAL(J)*PHI1+(1.-F_CAL(J))*PHI2)
C
      ENDDO
C
      GF1=SUM1
      GF2=MAX_1
      GF3=MAX_1-MIN_1
      GF4=SUM4/X_POINTS
      GF5=SUM5/X_POINTS
      GF6=SUM6/X_POINTS
      GF7=SUM7/X_POINTS
      GF8=SUM8/(ALPHA*(1.+ALPHA))
      GF9=SUM9
      GF10=SUM10/(X_POINTS*ALPHA*(ALPHA-1.0))
      GF11=SUM11/(ALPHA*(ALPHA+1.0))
      GF12=SUM12
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE KERNEL(EXPE,CALC,W,N_POINTS,ALPHA,C,D,CNORM,
     1                  KD1,KD2,KD3,KD4,KD5,KD6,KD7,KD8,KD9,KD10,
     2                  KD11,KD12)
C
C   This subroutine computes a kernel distance between two data sets
C
C     KD1  :  Linear kernel
C     KD2  :  Polynomial kernel
C     KD3  :  Sigmoid kernel
C     KD4  :  Hypertangent kernel
C     KD5  :  Exponential kernel
C     KD6  :  Gaussian kernel
C     KD7  :  Rational quadratic kernel
C     KD8  :  Multiquadric kernel
C     KD9  :  Inverse multiquadric kernel
C     KD10 :  Cauchy kernel
C     KD11 :  Generalized T-student kernel
C     KD12 :  Linear spline kernel
C
C  Input parameters:
C
C       EXPE     :  y coordinates of the experimental file
C       CALC     :  y coordinates of the calculation file
C       W        :  weight function
C       N_POINTS :  number of points in the files
C       ALPHA    :  alpha parameter 
C       C        :  c/sigma parameter
C       D        :  d parameter
C       CNORM    :  scaling coefficient to rescale calculation to experiment
C
C
C  Output parameters:
C
C       KDn      :  kernel distance value
C
C
C   Author :  D. Sébilleau
C
C                                          Last modified : 29 Aug 2014
C
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 EXPE(N_SIZE),CALC(N_SIZE),W(N_SIZE)
      REAL*4 KD1,KD2,KD3,KD4,KD5,KD6,KD7,KD8,KD9,KD10,KD11,KD12
      REAL*4 ALPHA,C,SIGMA2
      REAL*4 SCALAR_11,SCALAR_22,SCALAR_12
      REAL*4 NORM2_11,NORM2_22,NORM2_12
      REAL*4 K1_11,K1_22,K1_12
      REAL*4 K2_11,K2_22,K2_12
      REAL*4 K3_11,K3_22,K3_12
      REAL*4 K4_11,K4_22,K4_12
      REAL*4 K5_11,K5_22,K5_12
      REAL*4 K6_11,K6_22,K6_12 
      REAL*4 K7_11,K7_22,K7_12 
      REAL*4 K8_11,K8_22,K8_12 
      REAL*4 K9_11,K9_22,K9_12 
      REAL*4 K10_11,K10_22,K10_12 
      REAL*4 K11_11,K11_22,K11_12
      REAL*4 K12_11,K12_22,K12_12
      REAL*4 CNORM,SCALING
      REAL*4 MINXY(N_SIZE),X,Y,Z,SCAL
C
      INTEGER D
C
      SIGMA2=C*C
C
      SCALING=1./FLOAT(N_POINTS*N_POINTS)
      SCAL=1.5*SQRT(FLOAT(N_POINTS))
C
      SCALAR_11=0.
      SCALAR_22=0.
      SCALAR_12=0.
C
      NORM2_11=0.
      NORM2_22=0.
      NORM2_12=0.
C
C  Computing <X|Y> and ||X-Y||^2
C
      DO I=1,N_POINTS
C
        SCALAR_11=SCALAR_11+EXPE(I)*EXPE(I)
        SCALAR_22=SCALAR_22+CNORM*CNORM*CALC(I)*CALC(I)
        SCALAR_12=SCALAR_12+EXPE(I)*CNORM*CALC(I)
C
        NORM2_12=NORM2_12+(EXPE(I)-CNORM*CALC(I))*
     1                    (EXPE(I)-CNORM*CALC(I))
C
      ENDDO
C
C  Scaling by (1/N_POINTS^2) to prevent tanh, exp(- ...) , etc 
C       from having a very large argument which will prevent
C       any difference to appear between experiment and calculation
C
      SCALAR_11=SCALAR_11*SCALING
      SCALAR_22=SCALAR_22*SCALING
      SCALAR_12=SCALAR_12*SCALING
      NORM2_12=NORM2_12*SCALING
C
C  Computing the minimum of EXPE and CALC for spline kernel
C
      DO I=1,N_POINTS
C
        MINXY(I)=MIN(EXPE(I),CNORM*CALC(I))
C
      ENDDO
C
C  Linear kernel
C
      K1_11=SCALAR_11
      K1_22=SCALAR_22
      K1_12=SCALAR_12
C
C  Polynomial kernel
C
      K2_11=(ALPHA*SCALAR_11+C)**D
      K2_22=(ALPHA*SCALAR_22+C)**D
      K2_12=(ALPHA*SCALAR_11+C)**D
C
C  Sigmoid kernel
C
      K3_11=TANH(ALPHA*SCALAR_11+C)
      K3_22=TANH(ALPHA*SCALAR_22+C)
      K3_12=TANH(ALPHA*SCALAR_11+C)
C
C  Hypertangent kernel
C
      K4_11=1.0-TANH(NORM2_11/SIGMA2)
      K4_22=1.0-TANH(NORM2_22/SIGMA2)
      K4_12=1.0-TANH(NORM2_12/SIGMA2)
C
C  Exponential kernel
C
      K5_11=EXP(-SQRT(NORM2_11)/(2.0*SIGMA2))
      K5_22=EXP(-SQRT(NORM2_22)/(2.0*SIGMA2))
      K5_12=EXP(-SQRT(NORM2_12)/(2.0*SIGMA2))
C
C  Gaussian kernel
C
      K6_11=EXP(-NORM2_11/(2.0*SIGMA2))
      K6_22=EXP(-NORM2_22/(2.0*SIGMA2))
      K6_12=EXP(-NORM2_12/(2.0*SIGMA2))
C
C  Rational quadratic kernel
C
      K7_11=1.0-(NORM2_11/(NORM2_11+SIGMA2))
      K7_22=1.0-(NORM2_22/(NORM2_22+SIGMA2))
      K7_12=1.0-(NORM2_12/(NORM2_12+SIGMA2))
C
C  Multiquadric kernel
C
      K8_11=SQRT(SIGMA2-NORM2_11)
      K8_22=SQRT(SIGMA2-NORM2_22)
      K8_12=SQRT(SIGMA2-NORM2_12)
C
C  Inverse multiquadric kernel
C
      K9_11=1.0/SQRT(NORM2_11+SIGMA2)
      K9_22=1.0/SQRT(NORM2_22+SIGMA2)
      K9_12=1.0/SQRT(NORM2_12+SIGMA2)
C
C  Cauchy kernel
C
      K10_11=1.0/(1.0+(NORM2_11/SIGMA2))
      K10_22=1.0/(1.0+(NORM2_22/SIGMA2))
      K10_12=1.0/(1.0+(NORM2_12/SIGMA2))
C
C  Generalized T-student kernel
C
      K11_11=1.0/(1.0+NORM2_11**D)
      K11_22=1.0/(1.0+NORM2_22**D)
      K11_12=1.0/(1.0+NORM2_12**D)
C
C  Linear spline kernel
C
      K12_11=1.0
      K12_22=1.0
      K12_12=1.0
C
      DO I=1,N_POINTS
C
        X=EXPE(I)
        Y=CALC(I)
        Z=MINXY(I)
        K12_11=K12_11*(1.0+X*X+X*X*X-0.5*(X+X)*X*X+X*X*X/3.)/SCAL
        K12_22=K12_22*(1.0+Y*Y+Y*Y*Y-0.5*(Y+Y)*Y*Y+Y*Y*Y/3.)/SCAL
        K12_12=K12_12*(1.0+X*Y+X*Y*Z-0.5*(X+Y)*Z*Z+Z*Z*Z/3.)/SCAL
C
      ENDDO      
C
      KD1=SQRT(K1_11+K1_22-2.*K1_12)
      KD2=SQRT(K2_11+K2_22-2.*K2_12)
      KD3=SQRT(K3_11+K3_22-2.*K3_12)
      KD4=SQRT(K4_11+K4_22-2.*K4_12)
      KD5=SQRT(K5_11+K5_22-2.*K5_12)
      KD6=SQRT(K6_11+K6_22-2.*K6_12)
      KD7=SQRT(K7_11+K7_22-2.*K7_12)
      KD8=SQRT(K8_11+K8_22-2.*K8_12)
      KD9=SQRT(K9_11+K9_22-2.*K9_12)
      KD10=SQRT(K10_11+K10_22-2.*K10_12)
      KD11=SQRT(K11_11+K11_22-2.*K11_12)
      KD12=SQRT(K12_11+K12_22-2.*K12_12)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE MOMENTS(I2,X2,N_POINTS,N_GRID,IUO1,I_ALG,MU,NU,BASIS,
     1                   M,N_MOM)
C
C  This subroutine computes the 1D moments of function I(X)
C      up to order N_MOM and stores them in array M. 
C      BASIS defines the type of moment computed:
C
C                  BASIS = GEOM  geometric moments
C                  BASIS = LEGE  continuous Legendre 
C                  BASIS = CHEB  discrete Chebyshev
C                  BASIS = KRAW  discrete Krawtchouk
C                  BASIS = HAHN  discrete Hahn
C                  BASIS = MEIX  discrete Meixner
C                  BASIS = CHAR  discrete Charlier
C                  BASIS = SHMA  discrete Shmaliy
C
C
C  Input parameters:
C
C       I2       :  y coordinates of the input file
C       X2       :  x coordinates of the input file
C       N_POINTS :  number of points in the file
C       N_GRID   :  number of grid points on which the basis functions 
C                     are defined 
C       IUO1     :  checkfile index for printing
C       I_ALG    :  switch to select the type of recurrence used
C       MU       :  mu parameter
C       NU       :  nu parameter
C       BASIS    :  type of basis functions used to compute the moments
C
C
C  Output parameters:
C
C       M        :  moment values
C       N_MOM    :  number of moments
C
C   Author :  D. Sébilleau
C
C                                          Last modified :  9 Sep 2014
C
      PARAMETER (N_SIZE=1000,NMAX=9999)
C
      REAL*4 I2(N_SIZE),X2(N_SIZE),M(0:NMAX)
      REAL*4 I3(N_SIZE),X3(N_SIZE)
      REAL*4 YS(N_SIZE)
C
      REAL*8 MU,NU
C
      CHARACTER*4 BASIS
C
      IF(BASIS.EQ.'GEOM') THEN
C
C     Continuous geometric moments
C
        CALL M_GEOMETRIC(I2,X2,N_POINTS,M,N_MOM)
C
      ELSEIF(BASIS.EQ.'LEGE') THEN
C
C     Continuous Legendre moments
C
        CALL M_LEGENDRE(I2,X2,N_POINTS,M,N_MOM)
C
      ELSE
C
C     Discrete polynomials moments on a uniform lattice grid
C
C
C  Step for interpolation on the lattice grid
C
        STEP_D=(X2(N_POINTS)-X2(1))/FLOAT(N_GRID-1)
C
C  Interpolation of the input curve
C
        DO J=1,N_GRID
          X3(J)=X2(1)+FLOAT(J-1)*STEP_D
        ENDDO     
C
        YP1=2.0E30
        YPN=2.0E30
C
        CALL SPLINE(X2,I2,N_POINTS,YP1,YPN,YS)
C
        DO J=1,N_GRID
          XX=X3(J)
          CALL SPLINT(X2,I2,YS,N_POINTS,XX,YY,*10)
          I3(J)=YY
        ENDDO
C
C
C     Discrete Chebyshev moments  
C
        IF(BASIS.EQ.'CHEB') THEN
C
          CALL M_CHEBYCHEV(I3,X3,N_GRID,M,N_MOM+1,*20)
C
C     Discrete Krawtchouk, Hahn, Meixner, Charlier or Shmaliy moments
C
        ELSE
C
          CALL M_ORTH_POLY(BASIS,I_ALG,I3,X3,MU,NU,N_GRID,
     1                     M,N_MOM,*20)
        ENDIF 
C
      ENDIF
      GOTO 5
C
  10  WRITE(IUO1,11)
      STOP
  20  WRITE(IUO1,21) NMAX
      STOP
C
C  Formats
C
  11  FORMAT(//,10X,'<<<<<  WRONG VALUE FOR XA IN SUBROUTINE ',
     1       '   SPLINT  >>>>>',//)                                                  
  21  FORMAT(//,10X,'-->  DATA ERROR : N_MOM SHOULD BE ',
     1      'LOWER THAN ',I9,//)
C
   5  RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE M_GEOMETRIC(I,X,N_POINTS,M,N_MOM)
C
C  This subroutine computes the 1D geometric moments of function I(X)
C                        up to order N_MOM,
C
C
C  Input parameters:
C
C       I        :  y coordinates of the input file
C       X        :  x coordinates of the input file
C       N_POINTS :  number of points in the file
C
C
C  Output parameters:
C
C       M        :  moment values
C       N_MOM    :  number of moments
C
C
C   Author :  D. Sébilleau
C
C                                          Last modified : 30 Jul 2014
C
      PARAMETER (N_SIZE=1000,NMAX=9999)
C
      REAL*4 I(N_SIZE),X(N_SIZE),M(0:NMAX)
C
      REAL*8 SUM_X,DELTA,A,B
C
      INTEGER P
C
      DELTA=DBLE(X(2)-X(1))
C
      DO P=0,N_MOM
        
C
        SUM_X=0.D0
        DO J=1,N_POINTS
C
          A=(DBLE(X(J))+DELTA*0.5D0)**(P+1)
          B=(DBLE(X(J))-DELTA*0.5D0)**(P+1)
          SUM_X=SUM_X+DBLE(I(J))*(A-B)/DFLOAT(P+1) 
C
        ENDDO
C
        M(P)=REAL(SUM_X)
        WRITE(98,*) P,M(P)
        
C
      ENDDO
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE M_LEGENDRE(I,X,N_POINTS,MGL,N_MOM)
C
C  This subroutine computes the 1D Legendre moments MGL of function I(X)
C       up to order N_MOM, and writes the reconstructed function 
C        I_REC for checking. X is an angle expressed in degrees
C
C  The integration is performed using a NP-points Gauss-Legendre method
C
C
C
C  Input parameters:
C
C       I        :  y coordinates of the input file
C       X        :  x coordinates of the input file
C       N_POINTS :  number of points in the file
C
C
C  Output parameters:
C
C       MGL      :  moment values
C       N_MOM    :  number of moments
C
C
C   Authors :  D. Sébilleau and K. Hatada  
C
C  First version: October 2012             Last modified : 20 Jun 2014
C
      PARAMETER (N_SIZE=1000,NMAX=9999,NP=250)
C
      REAL*4 I(N_SIZE),X(N_SIZE),PL(0:NMAX)
      REAL*4 SUM_X,SUM_L
      REAL*4 ALM(0:NMAX,0:NMAX),SUM_K
      REAL*4 XGL(NP),WGT(NP),X1,X2,FUNC(NP),YS(N_SIZE),XRAD(N_SIZE)
      REAL*4 XX(NP),MGL(0:NMAX),MGGL(0:NMAX),LEG(0:NMAX,NP)
      REAL*4 HFBICO
C
      INTEGER P
C
      DATA PI /3.141593/
C
C Mapping of the original X points onto an interval more suitable 
C    for the Legendre integration (i.e. avoiding the range around 
C    zero degrees where there are unstabilities)
C
      SFT=50.
      SCL=0.5
      DO J = 1, N_POINTS
        X(J)=(X(J)+SFT)*SCL
      ENDDO
C
C Changing angles to radians
C  
      DO J = 1, N_POINTS
        XRAD(J)=X(J)*PI/180.
      ENDDO
C
C  X1 is the minimum of region of integral, and x2 is the maximum of integral
C
      X1=XRAD(1)
      X2=XRAD(N_POINTS)
C
C  Construct Gauss-Legendre points from Numerical Recipes subroutine
C
      CALL GAULEG(X1,X2,XGL,WGT,NP)
C
C  Natural Spline interpolation of f(x) at the Gauss-Legendre points
C
      YP1=2.0E20
      YPN=2.0E20
C
      CALL SPLINE(XRAD,I,N_POINTS,YP1,YPN,YS)
C
      DO K = 1, NP
        XX(K)=XGL(K)
        CALL SPLINT(XRAD,I,YS,N_POINTS,XX(K),YY,*6)
        FUNC(K)=YY
      ENDDO
C
      GOTO 9
C
  6   WRITE(66,7)
C
  9   CONTINUE
C
C  Calculation of the coefficient (Legendre polynomial Wikipedia JP)
C
      ALM(0,0)=1.0
      DO J=1,N_MOM
        ALM(0,J)=0.0
      ENDDO
      DO K=1,N_MOM
        DO J=0,K
          IF (MOD(K+J-1,2).EQ.0.AND.K+J-1.NE.0.AND.
     1        (K+J-1)/2.GE.K) THEN
            ALM(K,J)=REAL(2**K)*BICO(K,J)*BICO((K+J-1)/2,K)
          ELSE IF (MOD(K+J-1,2).NE.0.AND.K+J-1.NE.0) THEN
            ALM(K,J)=REAL(2**K)*BICO(K,J)*HFBICO(K+J-1,K)
          ELSE
            ALM(K,J)=0.0
          ENDIF
        ENDDO
        DO J=K+1,N_MOM
          ALM(K,J)=0.0
        ENDDO
      ENDDO
C
C  Storage of the Legendre polynomials for Gauss-Legendre
C
      IF (.TRUE.) THEN
C
C  Construction of Legendre function by recursion
C
        DO J=1,NP
C
          ANGLE=XX(J)
          CALL POLLEG(N_MOM+1,COS(ANGLE),PL)
C
          DO L=0,N_MOM
            LEG(L,J)=PL(L)
          ENDDO
C
        ENDDO

      ELSE
C
C  Construction of Legendre function by power series
C
        DO J=1,NP
C
          ANGLE=XX(J)
C
          DO L=0,N_MOM
            SUM_X=0.0
            IF (ABS(COS(XX(J))).LE.1.0) THEN
              DO K=L,0,-1
                SUM_X=SUM_X+ALM(L,K)*COS(XX(J))**K
              ENDDO
            ELSE
              DO K=0,L
                SUM_X=SUM_X+ALM(L,K)*COS(XX(J))**K
              ENDDO
            END IF
            LEG(L,J)=SUM_X
          ENDDO
C
        ENDDO  
C 
      ENDIF
C
C  Computation of the moments for Gauss-Legendre
C
      DO L=0,N_MOM
C
        SUM_X=0.
        DO J=1,NP
C
          SUM_X=SUM_X+FUNC(J)*LEG(L,J)*SIN(XX(J))*WGT(J)
C
        ENDDO
C
        MGL(L)=FLOAT(L+L+1)*0.5*SUM_X
        WRITE(98,*) L,MGL(L)
C
      ENDDO
C
C  Computation of the geometric moments by Gauss-Legendre integration
C
      DO P=0,N_MOM
C
        SUM_X=0.
        DO J=1,NP
C
          SUM_X=SUM_X+FUNC(J)*COS(XX(J))**P*SIN(XX(J))*WGT(J)
C
        ENDDO
C
        MGGL(P)=SUM_X
C
      ENDDO
C
C  Check: computing Legendre moments from geometric ones 
C
      DO L=0,N_MOM
        SUM_K=0.
        DO K=0,L
          SUM_K=SUM_K+ALM(L,K)*MGGL(K)
        ENDDO
        SUM_K=SUM_K*FLOAT(L+L+1)*0.5
C

      ENDDO

C
C  Reconstruction of the original function:
C
C     SUM_L : Legendre moments computed directly (MGL)
C     SUM_LM: Legendre moments computed from geometric moments (SUM_K)
C
      DO J=1,NP
C
        SUM_L=0.
        SUM_LM=0.
        DO L=0,N_MOM
C
          SUM_L=SUM_L+MGL(L)*LEG(L,J)
          SUM_K=0.
          DO K=0,L
            SUM_K=SUM_K+ALM(L,K)*MGGL(K)
          ENDDO
          SUM_K=SUM_K*FLOAT(2*L+1)*0.5
          SUM_LM=SUM_LM+SUM_K*LEG(L,J)
C
        ENDDO
C
C  Transforms back the X values to the original interval
C
        WRITE(99,*) (XX(J)/PI*180/SCL-SFT),SUM_L,SUM_LM
C
      ENDDO
C
   7  FORMAT(//,10X,'<<<<<  WRONG VALUE FOR XA IN SUBROUTINE ',
     1       '   SPLINT  >>>>>',//)                                             
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE M_CHEBYCHEV(I,X,N_POINTS,MCH2,N_MAX,*)
C
C  This subroutine computes the 1D discrete Chebyshev moments MCH 
C       of function I(X) up to order N_MAX, and writes the 
C       reconstructed function I_REC for checking. X is an angle 
C       expressed in degrees
C
C   Based on: G. Wang and S. Wang, Pattern Recognition 39, 47-56 (2006)
C
C
C  Input parameters:
C
C       I        :  y coordinates of the input file
C       X        :  x coordinates of the input file
C       N_POINTS :  number of points in the file
C
C
C  Output parameters:
C
C       MCH2     :  moment values
C       N_MOM    :  number of moments
C
C
C   Author :  D. Sébilleau 
C
C                                          Last modified : 14 Jun 2013
C
      PARAMETER (N_SIZE=1000,NMAX=9999)
C
      REAL*4 I(N_SIZE),X(N_SIZE)
      REAL*4 MCH2(0:NMAX)
C
      REAL*8 I2(0:N_SIZE)
      REAL*8 TN(0:NMAX),RHO_N(0:NMAX),MCH(0:NMAX)
      REAL*8 ALPHA_A(0:NMAX,0:N_POINTS-1),BETA_A(0:N_POINTS)
      REAL*8 ALPHA_B(0:NMAX,0:N_POINTS),BETA_B(0:N_POINTS)
      REAL*8 PSI_A(0:NMAX,0:N_POINTS+1)
      REAL*8 PSI_B(0:NMAX+1,0:N_POINTS-1)
      REAL*8 X_POINTS,XX,XN,PROD_N,DENOM
C
      X_POINTS=DFLOAT(N_POINTS)
C
C  Checking the consistency: N_MAX <= N_POINTS
C    Otherwise, the normalization coefficient RHO_N is zero
C
      IF(N_MAX.GE.N_POINTS) THEN
        N_MAX=N_POINTS
        WRITE(6,10) N_POINTS
      ENDIF
C
C  Shifting I by -1
C
      DO J=1,N_POINTS
        I2(J-1)=DBLE(I(J))
      ENDDO
C
C  Computing the scaled Chebyshev polynomials at J=0: TN(N): eq (13)
C
      TN(0)=1.D0
      DO N=1,N_MAX-1
        XN=DFLOAT(N)
        TN(N)=TN(N-1)*(XN-X_POINTS)/X_POINTS
      ENDDO      
C
C  Computing the alpha_a and beta_a coefficients used for the moments: eq (15)
C     !!!!! Warning: misprint for alpha: 2(x+1) should be (2x+1) !!!!!
C
      DO J=0,N_POINTS-2
        XX=DFLOAT(J)
        DENOM=((XX+1.D0)*(X_POINTS-XX-1.D0))
        BETA_A(J)=XX*(XX-X_POINTS)/DENOM
        DO N=0,N_MAX-1
          XN=DFLOAT(N)
          ALPHA_A(N,J)=-(XN*(XN+1.D0)+(XX+XX+1.D0)*(XX-X_POINTS)+
     1                   XX+1.D0)/DENOM
        ENDDO
      ENDDO
C
C                             eq (37)
C
      BETA_A(N_POINTS-1)=0.D0
      DO N=0,N_MAX-1
        ALPHA_A(N,N_POINTS-1)=0.D0
      ENDDO
      BETA_A(N_POINTS)=0.D0
C
C  Computing the alpha_b and beta_b coefficients used to reconstruct:
C
      DO N=0,N_MAX-1
        XN=DFLOAT(N)
        DENOM=(XN+1.D0)*X_POINTS
        BETA_B(N)=-XN*(X_POINTS*X_POINTS-XN*XN)/
     1                (DENOM*X_POINTS)
        DO J=0,N_POINTS-1
          XX=DFLOAT(J)
          ALPHA_B(N,J)=(XN+XN+1.D0)*(XX+XX-X_POINTS+1.D0)/DENOM
        ENDDO
      ENDDO
C
      BETA_B(N_MAX)=0.D0
C
C  Normalization coefficient RHO_N: eq (3)
C
      PROD_N=X_POINTS
      RHO_N(0)=X_POINTS
      DO N=1,N_MAX-1
        XN=DFLOAT(N)
        PROD_N=PROD_N*(1.D0-(XN/X_POINTS)*(XN/X_POINTS))
        RHO_N(N)=PROD_N/(XN+XN+1.D0)
      ENDDO
C
C  Psi_a functions --> we need only PSI_A(N,0): eq (21)
C
      DO N=0,N_MAX-1
        XN=DFLOAT(N)
        PSI_A(N,N_POINTS+1)=0.D0
        PSI_A(N,N_POINTS)=0.D0
        DO J=N_POINTS-1,0,-1
          PSI_A(N,J)=ALPHA_A(N,J)*PSI_A(N,J+1)+BETA_A(J+1)*
     1                            PSI_A(N,J+2)+I2(J)
        ENDDO
      ENDDO
C
C  Chebychev moment
C
      DO N=0,N_MAX-1
        MCH(N)= PSI_A(N,0)*TN(N)/RHO_N(N)
        MCH2(N)=REAL(MCH(N))
        WRITE(98,*) N,MCH2(N)
      ENDDO
C
C  Psi_b functions --> we need only PSI_B(0,J): eq (39)
C
      DO J=0,N_POINTS-1
        PSI_B(N_MAX+1,J)=0.D0
        PSI_B(N_MAX,J)=0.D0
        DO N=N_MAX-1,0,-1
          PSI_B(N,J)=ALPHA_B(N,J)*PSI_B(N+1,J)+BETA_B(N+1)*PSI_B(N+2,J)
     1                                        +MCH(N)
        ENDDO
      ENDDO
C
C  Reconstruction of the original function: PSI_B(0,J): eq (41)
C
      DO J=0,N_POINTS-1
        WRITE(99,*) X(J+1),REAL(PSI_B(0,J))
      ENDDO
C
      RETURN
C
   10 FORMAT(5X,'!!!!!  N_MAX EXCEEDS THE NUMBER OF POINTS N_POINTS = '
     1       ,I5,' OF THE CURVE  !!!!!',/,5X,
     2        '!!!!!                   N_MAX TRUNCATED TO N_POINTS-1',
     3        '                    !!!!!',//)
C
      END
C
C=======================================================================
C
      SUBROUTINE M_ORTH_POLY(ORTH_POL,I_ALG,I,X,MU,NU,N_POINTS,
     1                       MCH2,N_MAX,*)
C
C  This subroutine computes the 1D discrete orthogonal polynomial
C       moments MCH of function I(X) up to order N_MAX, and writes 
C       the reconstructed function I_REC for checking. It assumes 
C       that I(X) is given on a uniform lattice grid
C
C  It considers Chebyshev, Krawtchouk, Hahn, Meixner, Charlier 
C       and Shmaliy polynomials. It is based on:
C
C     H. Zhu, M. Liu, H. Shu, H. Zhang and L. Luo,
C        IET Image Process., 2010, Vol. 4, Iss. 5, pp. 335–352
C
C     P. Ananth Raj and  A. Venkataramana, 
C        Image Processing, 2007. ICIP 2007. 
C        IEEE International Conference on Image Processing, pp. 37-40 
C        for Krawtchouk polynomials
C
C     L. J. Morales-Mendoza, H. Gamboa-Rosales and Y. S. Shmaliy,
C        Signal Processing 93, 1785-1793 (2013)
C
C  A summary of these discrete orthogonal polynomials is available 
C        in the MsSpec working notes
C
C  Input/output parameters :
C
C     ORTH_POL='CHEB' : Chebychev polynomials
C     ORTH_POL='KRAW' : Krawtchouk polynomials
C     ORTH_POL='HAHN' : Hahn polynomials
C     ORTH_POL='MEIX' : Meixner polynomials
C     ORTH_POL='CHAR' : Charlier polynomials
C     ORTH_POL='SHMA' : Shmaliy polynomials
C
C     I_ALG=1         : n-recurrence scheme
C     I_ALG=2         : x-recurrence scheme
C     I_ALG=3         : Clenshaw recurrence scheme
C
C     I(X)            : value of the input function at point X
C
C     MU, NU          : additional parameters used to define the 
C                       discrete orthogonal polynomials (see MsSpec 
C                       working notes)
C                           CHEB : not used
C                           KRAW : MU = p
C                           HAHN : MU = a and NU = b
C                           MEIX : MU = beta and NU = mu
C                           CHAR : MU = a_1
C                           SHMA : not used
C
C     MCH2            : moments of I(X)
C
C     N_POINTS        : number of X points of I(X)
C     N_MAX           : number of moments computed
C                       The moments vary from 0 to N_MAX
C
C    WARNING: the input I(X) has X varying from 1 to N_POINTS while 
C             the discrete orthogonal polynomials are defined on a 
C             uniform grid varying from 0 to (N_POINTS-1). Therefore,
C             a systematic shift is done in the formulas given by 
C             the previous articles or the MsSpec working notes
C
C   Author:  D. Sébilleau 
C
C                                          Last modified : 2 May 2014
C
      PARAMETER (N_SIZE=1000,NMAX=9999)
C
      REAL*4 I(N_SIZE),X(N_SIZE)
      REAL*4 MCH2(0:NMAX)
C
      REAL*8 ON(0:NMAX,0:N_POINTS-1),MCH(0:NMAX)
      REAL*8 A(0:NMAX),B(0:NMAX),C(0:NMAX),D(0:NMAX),E(0:NMAX)
      REAL*8 SIG(0:N_POINTS-1),TAU(0:N_POINTS-1),W(0:N_POINTS-1)
      REAL*8 X_POINTS,XX,SUM_X,MU,NU,CF1,CF2,CF3
C
      CHARACTER*4 ORTH_POL
C
      X_POINTS=DFLOAT(N_POINTS)
C
C  Checking the consistency: N_MAX <= N_POINTS
C    (except for Meixner and Charlier polynomials)
C
      IF((ORTH_POL.NE.'MEIX').AND.(ORTH_POL.NE.'CHAR')) THEN
        IF(N_MAX.GE.N_POINTS) THEN
          N_MAX=N_POINTS
          WRITE(6,10) N_POINTS
        ENDIF
      ENDIF
C
C  Shmaliy recurrence only defined in n
C
      IF(ORTH_POL.EQ.'SHMA') THEN
        I_ALG=1
      ENDIF
C
      IF(I_ALG.EQ.1) THEN
C
C  n-recurrence algorithm
C
C
C  Computing the n-recurrence coefficients: Table 2
C
        CALL NREC_COEF(X_POINTS,MU,NU,N_MAX,ORTH_POL,A,B,C,D,E)
C
C  Computing the initial values of the n-recurrence
C
        CALL NREC_INIT(MU,NU,N_POINTS,ORTH_POL,ON)
C
C  Computing the scaled discrete orthogonal polynomials OP(N,J) 
C                by the n-recursion: eq (9)
C
C  Coefficient B(N) here is (XX-B) where B is the coefficient 
C                occuring in eq (9)
C
        DO J=0,N_POINTS-1
          XX=DBLE(J)
          DO N=2,N_MAX
C
            IF(ORTH_POL.EQ.'KRAW') THEN
              ON(N,J)=A(N)*(B(N)-XX)*ON(N-1,J)-C(N)*ON(N-2,J)
            ELSEIF(ORTH_POL.EQ.'SHMA') THEN
              ON(N,J)=(A(N)+B(N)*XX)*D(N)*ON(N-1,J)+C(N)*E(N)*ON(N-2,J)
            ELSE
              ON(N,J)=((XX-B(N))*D(N)*ON(N-1,J)+C(N)*E(N)*ON(N-2,J))/
     1                A(N)
            ENDIF
C
          ENDDO
        ENDDO
C
      ELSEIF(I_ALG.EQ.2) THEN
C
C  x-recurrence algorithm
C
C
C  Computing the x-recurrence coefficients: Table 1
C
        CALL XREC_COEF(N_POINTS,MU,NU,N_MAX,ORTH_POL,SIG,TAU,W,A,B)
C
C  Computing the initial values of the x-recurrence
C
        CALL XREC_INIT(N_POINTS,MU,NU,N_MAX,ORTH_POL,W,B,ON)
C
C  Computing the scaled discrete orthogonal polynomials OP(N,J) 
C                by the x-recursion: eq (27)
C
        DO N=0,N_MAX
C
          DO J=2,N_POINTS-1
C
            CF1=(2.D0*SIG(J-1)+TAU(J-1)-A(N))/SQRT(W(J-1))
            CF2=-SIG(J-1)/SQRT(W(J-2))
            CF3=(SIG(J-1)+TAU(J-1))/SQRT(W(J))
            ON(N,J)=(CF1*ON(N,J-1)+CF2*ON(N,J-2))/CF3
          ENDDO
        ENDDO
C
      ENDIF
C
C  Discrete orthogonal polynomial moment: eq (28)
C
      DO N=0,N_MAX
        SUM_X=0.D0
        DO J=0,N_POINTS-1
          SUM_X=SUM_X+ON(N,J)*DBLE(I(J+1))
        ENDDO
        MCH(N)=SUM_X
        MCH2(N)=REAL(MCH(N))
        WRITE(98,*) N,MCH2(N)
      ENDDO
C
C  Reconstruction of the original function: eq (29)
C
      DO J=0,N_POINTS-1
        SUM_N=0.
        DO N=0,N_MAX-1
          SUM_N=SUM_N+REAL(MCH(N)*ON(N,J))
        ENDDO
        WRITE(99,*) X(J+1),SUM_N
      ENDDO
C
      RETURN
C
   10 FORMAT(5X,'!!!!!  N_MAX EXCEEDS THE NUMBER OF POINTS N_POINTS = '
     1       ,I5,' OF THE CURVE  !!!!!',/,5X,
     2        '!!!!!                   N_MAX TRUNCATED TO N_POINTS-1',
     3        '                    !!!!!',//)
C
      END
C
C=======================================================================
C
      SUBROUTINE NREC_COEF(X_POINTS,MU,NU,N_MAX,ORTH_POL,A,B,C,D,E)
C
C  This subroutine computes the n-recurrence coefficients for 
C    the different discrete orthogonal coeffcients from Table 2 of
C
C     H. Zhu, M. Liu, H. Shu, H. Zhang and L. Luo,
C        IET Image Process., 2010, Vol. 4, Iss. 5, pp. 335–352
C
C     P. Ananth Raj and  A. Venkataramana, 
C        Image Processing, 2007. ICIP 2007. 
C        IEEE International Conference on Image Processing, pp. 37-40 
C        for Krawtchouk polynomials (because of misprints is the 
C        previous article
C
C     L. J. Morales-Mendoza, H. Gamboa-Rosales and Y. S. Shmaliy,
C        Signal Processing 93, 1785-1793 (2013)
C
C
C  Input parameters:
C
C       X_POINTS   :  number of points in the function file
C       MU, NU     :  additional parameters used to define the 
C                       discrete orthogonal polynomials (see MsSpec 
C                       working notes)
C                           CHEB : not used
C                           KRAW : MU = p
C                           HAHN : MU = a and NU = b
C                           MEIX : MU = beta and NU = mu
C                           CHAR : MU = a_1
C                           SHMA : not used
C       N_MAX      :  number of moments computed
C       ORTH_POL   :  type of orthogonal polynomials used
C
C                   ='CHEB' : Chebychev polynomials
C                   ='KRAW' : Krawtchouk polynomials
C                   ='HAHN' : Hahn polynomials
C                   ='MEIX' : Meixner polynomials
C                   ='CHAR' : Charlier polynomials
C                   ='SHMA' : Shmaliy polynomials
C
C
C  Output parameters:
C
C       A,B,C,D,E  :  coefficients for the n-recurrence 
C
C
C   Author:  D. Sébilleau 
C
C                                          Last modified :  2 May 2014
C
      PARAMETER (NMAX=9999)
C
      REAL*8 A(0:NMAX),B(0:NMAX),C(0:NMAX),D(0:NMAX),E(0:NMAX)
      REAL*8 X_POINTS,MU,NU,XN,XN1,XN2,XN3,NUM,DEN,COEF
C
      CHARACTER*4 ORTH_POL
C
      DO N=2,N_MAX
C
        XN=DFLOAT(N)
        XN1=XN-1.D0
        XN2=XN+XN-1.D0
        XN3=XN+XN+1.D0
C
        IF(ORTH_POL.EQ.'CHEB') THEN
C
C  Discrete Chebyshev polynomials
C
          A(N)=0.5D0*XN/XN2
          B(N)=(X_POINTS-1.D0)*0.5D0
          C(N)=-0.5D0*XN1*(X_POINTS*X_POINTS-XN1*XN1)/XN2
          D(N)=SQRT(XN3/((X_POINTS*X_POINTS-XN*XN)*XN2))
          E(N)=SQRT(XN3/((X_POINTS*X_POINTS-XN*XN)*
     1                   (X_POINTS*X_POINTS-XN1*XN1)*(XN2-2.D0)))
C
        ELSEIF(ORTH_POL.EQ.'KRAW') THEN
C
C  Discrete Krawtchouk polynomials (using the Raj-Venkataramana 
C          article)
C
          A(N)=DSQRT(1.D0/(MU*(1.D0-MU)*XN*(X_POINTS-XN)))
          B(N)=MU*(X_POINTS-2.D0*XN+1.D0)+XN1
          C(N)=DSQRT((XN1*(X_POINTS-XN1))/((X_POINTS-XN)*XN))
          D(N)=0.D0
          E(N)=0.D0
C
        ELSEIF(ORTH_POL.EQ.'HAHN') THEN
C
C  Discrete Hahn polynomials
C
          NUM=XN*(MU+NU+XN)
          DEN=(MU+NU+XN2)*(MU+NU+XN2+1.D0)
          A(N)=NUM/DEN
C
          NUM=(NU*NU-MU*MU)*(MU+NU+X_POINTS+X_POINTS)
          DEN=(MU+NU+XN2-1.D0)*(MU+NU+XN2+1.D0)
          B(N)=0.25D0*((MU-NU+X_POINTS+X_POINTS-2.D0)+NUM/DEN)
C
          NUM=(MU+XN1)*(NU+XN1)*(MU+NU+X_POINTS+XN1)*(X_POINTS-XN1)
          DEN=(MU+NU+XN2-1.D0)*(MU+NU+XN2)
          C(N)=-NUM/DEN
C
          NUM=XN*(MU+NU+XN)*(MU+NU+XN3)
          DEN=(X_POINTS-XN)*(MU+XN)*(NU+XN)*(MU+NU+XN2)*
     1        (MU+NU+XN+X_POINTS)
          D(N)=SQRT(NUM/DEN)
C
          NUM=XN*XN1*(MU+NU+XN)*(MU+NU+XN1)*(MU+NU+XN3)
          DEN=(MU+XN)*(MU+XN1)*(NU+XN)*(NU+XN1)*(X_POINTS-XN1)*
     1        (X_POINTS-XN)*(MU+NU+XN2-2.D0)*(MU+NU+XN+X_POINTS)*
     2        (MU+NU+XN1+X_POINTS)
          E(N)=SQRT(NUM/DEN)
C
        ELSEIF(ORTH_POL.EQ.'MEIX') THEN
C
C  Discrete Meixner polynomials
C
          A(N)=NU/(NU-1.D0)
          B(N)=(XN1+NU*XN-NU+MU*NU)/(1.D0-NU)
          C(N)=XN1*(XN-2.D0+MU)/(1.D0-NU)
          D(N)=SQRT(NU/(XN*(MU+XN1)))
          E(N)=SQRT(NU*NU/(XN*XN1*(MU+XN-2.D0)*(MU+XN1)))
C
        ELSEIF(ORTH_POL.EQ.'CHAR') THEN
C
C  Discrete Charlier polynomials
C
          A(N)=-MU
          B(N)=XN1+MU
          C(N)=XN1
          D(N)=SQRT(MU/XN)
          E(N)=SQRT(MU*MU/(XN*XN1))
C  
        ELSEIF(ORTH_POL.EQ.'SHMA') THEN
C
C  Discrete Shmaliy polynomials
C     (see L. J. Morales-Mendoza, H. Gamboa-Rosales 
C      and Y. S. Shmaliy
C
C                A(N) = alpha_n
C                B(N) = beta_n
C                C(N) = - zeta_n
C                D(N) = sqrt(d_{n-1}^2 / d_{n}^2)
C                E(N) = sqrt(d_{n-2}^2 / d_{n}^2)
C
C
          COEF=XN2*(X_POINTS+XN)
C
          A(N)=2.D0*XN*(2.D0*X_POINTS-1.D0)/COEF
          B(N)=-2.D0*(4.D0*XN-1.D0/XN)/COEF
          C(N)=-XN3*(X_POINTS-XN)/COEF
C
          NUM=XN*(X_POINTS-XN-2.D0)*(X_POINTS+XN)
          DEN=(XN+1.D0)*(X_POINTS-2.D0)*(X_POINTS-3.D0)
          D(N)=SQRT(NUM/DEN)
C
          NUM=XN1*(X_POINTS-XN-2.D0)*(X_POINTS-XN-3.D0)*(X_POINTS+XN)*
     1            (X_POINTS+XN1)
          DEN=(XN+1.D0)*(X_POINTS-2.D0)*(X_POINTS-3.D0)*(X_POINTS-4.D0)*
     1                  (X_POINTS-5.D0)
          E(N)=SQRT(NUM/DEN)
C
        ENDIF
C
      ENDDO
C
      END
C
C=======================================================================
C

      SUBROUTINE NREC_INIT(MU,NU,N_P,ORTH_POL,ON)
C
C  This subroutine computes the initial values for the n-recurrence
C     for the calculation of the discrete orthogonal polynomials
C
C     H. Zhu, M. Liu, H. Shu, H. Zhang and L. Luo,
C        IET Image Process., 2010, Vol. 4, Iss. 5, pp. 335–352
C
C     P. Ananth Raj and  A. Venkataramana, 
C        Image Processing, 2007. ICIP 2007. 
C        IEEE International Conference on Image Processing, pp. 37-40 
C        for Krawtchouk polynomials (because of misprints in the 
C        previous article) or equivalently
C     P.-T. Yap and R. Paramesran, IEEE Transactions on Image 
C        Processing, 12, 1367 (2003)
C
C     L. J. Morales-Mendoza, H. Gamboa-Rosales and Y. S. Shmaliy,
C        Signal Processing 93, 1785-1793 (2013)
C
C
C  Input parameters:
C
C       MU, NU     :  additional parameters used to define the 
C                       discrete orthogonal polynomials (see MsSpec 
C                       working notes)
C                           CHEB : not used
C                           KRAW : MU = p
C                           HAHN : MU = a and NU = b
C                           MEIX : MU = beta and NU = mu
C                           CHAR : MU = a_1
C                           SHMA : not used
C       N_P        :  number of points in the input function
C       ORTH_POL   :  type of orthogonal polynomials used
C
C                   ='CHEB' : Chebychev polynomials
C                   ='KRAW' : Krawtchouk polynomials
C                   ='HAHN' : Hahn polynomials
C                   ='MEIX' : Meixner polynomials
C                   ='CHAR' : Charlier polynomials
C                   ='SHMA' : Shmaliy polynomials
C
C
C  Output parameters:
C
C       ON         :  discrete orthogonal polynomial 
C
C
C   Author:  D. Sébilleau 
C
C                                          Last modified : 30 Apr 2014
C
C
      PARAMETER (NMAX=9999)
C
      REAL*8 ON(0:NMAX,0:N_P-1),XX,XP,XP1,MU,NU
      REAL*8 LNUM,LDEN,FACTLN,GAMMLN
C
      CHARACTER*4 ORTH_POL
C
      XP=DBLE(N_P)
      XP1=XP-1.D0
C
      DO J=0,N_P-1
C
        XX=DBLE(J)
C
        IF(ORTH_POL.EQ.'CHEB') THEN
C
C  Discrete Chebyshev polynomials
C
          ON(0,J)=SQRT(1.D0/XP)
          ON(1,J)=(XX+XX-XP1)*SQRT(3.D0/(XP*XP*XP-XP))
C
        ELSEIF(ORTH_POL.EQ.'KRAW') THEN
C
C  Discrete Krawtchouk polynomials
C     (see P. Ananth Raj and  A. Venkataramana or 
C       P.-T. Yap and R. Paramesran)
C
          LNUM=GAMMLN(XP)+XX*DLOG(MU)+(XP1-XX)*DLOG(1.D0-MU)
          LDEN=FACTLN(J)+GAMMLN(XP-XX)
          ON(0,J)=DEXP(0.5D0*(LNUM-LDEN))
C
          ON(1,J)=(1.D0-XX/(MU*XP1))*DSQRT(MU*XP1/(1.D0-MU))*ON(0,J)
C
        ELSEIF(ORTH_POL.EQ.'HAHN') THEN
C
C  Discrete Hahn polynomials
C
          LNUM=GAMMLN(XP+MU-XX)+GAMMLN(NU+XX+1.D0)+GAMMLN(XP)+
     1         GAMMLN(MU+NU+2.D0)
          LDEN=GAMMLN(XP-XX)+GAMMLN(XX+1.D0)+GAMMLN(MU+1.D0)+
     1         GAMMLN(NU+1.D0)+GAMMLN(MU+NU+XP+1.D0)
          ON(0,J)=DEXP(0.5D0*(LNUM-LDEN))
C
          LNUM=GAMMLN(XP+MU-XX)+GAMMLN(NU+XX+1.D0)+GAMMLN(XP1)+
     1         GAMMLN(MU+NU+2.D0)+DLOG(MU+NU+3.D0)
          LDEN=GAMMLN(XP-XX)+GAMMLN(XX+1.D0)+GAMMLN(MU+2.D0)+
     1         GAMMLN(NU+2.D0)+GAMMLN(MU+NU+XP+2.D0)
          ON(1,J)=((MU+NU+2.D0)*XX-(NU+1.D0)*XP1)*
     1             DEXP(0.5D0*(LNUM-LDEN))
C
        ELSEIF(ORTH_POL.EQ.'MEIX') THEN
C
C  Discrete Meixner polynomials
C
          LNUM=XX*DLOG(NU)+GAMMLN(MU+XX)+MU*DLOG(1.D0-NU)
          LDEN=FACTLN(J)+GAMMLN(MU)
          ON(0,J)=DEXP(0.5D0*(LNUM-LDEN))
C
          ON(1,J)=(MU+XX-XX/NU)*DSQRT(NU/MU)*ON(0,J)
C  
        ELSEIF(ORTH_POL.EQ.'CHAR') THEN
C
C  Discrete Charlier polynomials
C
          LNUM=-MU+XX*DLOG(MU)
          LDEN=FACTLN(J)
          ON(0,J)=DEXP(0.5D0*(LNUM-LDEN))
C
          ON(1,J)=(MU-XX)*ON(0,J)/DSQRT(MU)
C  
        ELSEIF(ORTH_POL.EQ.'SHMA') THEN
C
C  Discrete Shmaliy polynomials
C     (see L. J. Morales-Mendoza, H. Gamboa-Rosales 
C      and Y. S. Shmalyi
C
          ON(0,J)=DSQRT(2.D0*XX/(XP*XP1))
          ON(1,J)=DSQRT(XX*XP*(XP+1.D0)/(XP1*(XP-2.D0)))*
     1            (2.D0*(XP+XP1)-6.D0*XX)/(XP*(XP+1.D0))
C
        ENDIF
C
      ENDDO
C
      END
C
C=======================================================================
C
      SUBROUTINE XREC_COEF(N_P,MU,NU,N_MAX,ORTH_POL,SIG,TAU,W,
     1                     LAMBDA,DN2)
C
C  This subroutine computes the n-recurrence coefficients for 
C    the different discrete orthogonal coeffcients from Table 2 of
C
C     H. Zhu, M. Liu, H. Shu, H. Zhang and L. Luo,
C        IET Image Process., 2010, Vol. 4, Iss. 5, pp. 335–352
C
C     L. Zhu, J. Liao, X. Tong, L; Luo, B. Fu and G. Zhang,
C        Advances in Neural Networks, Lecture Notes in Computer Science,
C        Volume 5553, 310-317 (2009) (because of misprint in the 
C        previous article for Krawtchouk polynomials)
C
C
C  Input parameters:
C
C       N_P        :  number of points in the input function
C       MU, NU     :  additional parameters used to define the 
C                       discrete orthogonal polynomials (see MsSpec 
C                       working notes)
C                           CHEB : not used
C                           KRAW : MU = p
C                           HAHN : MU = a and NU = b
C                           MEIX : MU = beta and NU = mu
C                           CHAR : MU = a_1
C                           SHMA : not used
C       N_MAX      :  number of moments computed
C       ORTH_POL   :  type of orthogonal polynomials used
C
C                   ='CHEB' : Chebychev polynomials
C                   ='KRAW' : Krawtchouk polynomials
C                   ='HAHN' : Hahn polynomials
C                   ='MEIX' : Meixner polynomials
C                   ='CHAR' : Charlier polynomials
C                   ='SHMA' : Shmaliy polynomials
C
C
C  Output parameters:
C
C       SIG        :  sigma coefficient for x-recurrence 
C       TAU        :  tau coefficient for x-recurrence 
C       W          :  w coefficient for x-recurrence 
C       LAMBDA     :  lambda coefficient for x-recurrence 
C       DN2        :  dn square coefficient for x-recurrence 
C
C
C   Author:  D. Sébilleau 
C
C                                          Last modified :  2 May 2014
C
      REAL*8 SIG(0:N_P-1),TAU(0:N_P-1),W(0:N_P-1)
      REAL*8 LAMBDA(0:N_MAX),DN2(0:N_MAX),GAMMLN,FACTLN
      REAL*8 X_P,XX,XX1,MU,NU,XN,XN1,LNW,LDN2
C
      CHARACTER*4 ORTH_POL
C
      X_P=DBLE(N_P)
C
C  x-dependent coeffcients SIG, TAU and W
C
      DO J=0,N_P-1
C
        XX=DBLE(J)
        XX1=XX+1.D0
C
        IF(ORTH_POL.EQ.'CHEB') THEN
C
C  Discrete Chebyshev polynomials
C
          SIG(J)=XX*(X_P-XX)
          TAU(J)=X_P-1.D0-XX-XX
C
          W(J)=1.D0
C
        ELSEIF(ORTH_POL.EQ.'KRAW') THEN
C
C  Discrete Krawtchouk polynomials
C
          SIG(J)=XX*(1.D0-MU)
          TAU(J)=MU*X_P-XX
C
          LNW=GAMMLN(X_P+1.D0)+XX*DLOG(MU)+(X_P-XX)*DLOG(1.D0-MU)-
     1        GAMMLN(XX1)-GAMMLN(X_P-XX+1.D0)
          W(J)=DEXP(LNW)
C
        ELSEIF(ORTH_POL.EQ.'HAHN') THEN
C
C  Discrete Hahn polynomials
C
          SIG(J)=XX*(X_P+MU-XX)
          TAU(J)=(NU+1.D0)*(X_P-1.D0)-(MU+NU+2.D0)*XX
C
          LNW=GAMMLN(X_P+MU-XX)+GAMMLN(NU+XX1)-GAMMLN(X_P-XX)-
     1        GAMMLN(XX1)
          W(J)=DEXP(LNW)
C
        ELSEIF(ORTH_POL.EQ.'MEIX') THEN
C
C  Discrete Meixner polynomials
C
          SIG(J)=XX
          TAU(J)=MU*NU-XX*(1.D0-NU)
C
          LNW=XX*DLOG(NU)+GAMMLN(MU+XX)-GAMMLN(MU)-FACTLN(J)
          W(J)=DEXP(LNW)
C
        ELSEIF(ORTH_POL.EQ.'CHAR') THEN
C
C  Discrete Charlier polynomials
C
          SIG(J)=XX
          TAU(J)=MU-XX
C
          LNW=-MU+XX*DLOG(MU)-FACTLN(J)
          W(J)=DEXP(LNW)
C
        ENDIF
C
      ENDDO
C
C  n-dependent coefficients LAMBDA and DN2
C
      DO N=0,N_MAX-1
C
        XN=DBLE(N)
        XN1=XN+1.D0
C
        IF(ORTH_POL.EQ.'CHEB') THEN
C
          LAMBDA(N)=XN*(XN+1.D0)
C
          LDN2=FACTLN(N_P+N)-FACTLN(N_P-N-1)
          DN2(N)=DEXP(LDN2)/(XN+XN+1.D0)
C
        ELSEIF(ORTH_POL.EQ.'KRAW') THEN
C
          LAMBDA(N)=XN
C
          LDN2=XN*(DLOG(1.D0-MU)-DLOG(MU))+FACTLN(N)+FACTLN(N_P-N)-
     1         FACTLN(N_P)
          DN2(N)=DEXP(LDN2)
C
        ELSEIF(ORTH_POL.EQ.'HAHN') THEN

          LAMBDA(N)=XN*(MU+NU+XN1)

          LDN2=GAMMLN(MU+XN1)+GAMMLN(NU+XN1)+GAMMLN(MU+NU+XN1+X_P)-
     1         DLOG(MU+NU+XN+XN1)-FACTLN(N)-GAMMLN(X_P-XN)-
     2         GAMMLN(MU+NU+XN1)
          DN2(N)=DEXP(LDN2)

        ELSEIF(ORTH_POL.EQ.'MEIX') THEN
C
          LAMBDA(N)=XN*(1.D0-NU)
C
          LDN2=FACTLN(N)+GAMMLN(MU+N)-XN*DLOG(NU)-MU*DLOG(1.D0-NU)-
     1         GAMMLN(MU)
          DN2(N)=DEXP(LDN2)
C
        ELSEIF(ORTH_POL.EQ.'CHAR') THEN
C
          LAMBDA(N)=XN
C
          LDN2=FACTLN(N)-XN*DLOG(MU)
          DN2(N)=DEXP(LDN2)

        ENDIF
C
      ENDDO
C
      END
C
C=======================================================================
C

      SUBROUTINE XREC_INIT(N_P,MU,NU,N_MAX,ORTH_POL,W,DN2,ON)
C
C  This subroutine computes the initial values for the x-recurrence
C     for the calculation of the discrete orthogonal polynomials
C
C     H. Zhu, M. Liu, H. Shu, H. Zhang and L. Luo,
C        IET Image Process., 2010, Vol. 4, Iss. 5, pp. 335–352
C
C     L. Zhu, J. Liao, X. Tong, L; Luo, B. Fu and G. Zhang,
C        Advances in Neural Networks, Lecture Notes in Computer Science,
C        Volume 5553, 310-317 (2009) (because of misprint in the 
C        previous article for Krawtchouk polynomials)
C
C
C  Input parameters:
C
C       N_P        :  number of points in the input function
C       MU, NU     :  additional parameters used to define the 
C                       discrete orthogonal polynomials (see MsSpec 
C                       working notes)
C                           CHEB : not used
C                           KRAW : MU = p
C                           HAHN : MU = a and NU = b
C                           MEIX : MU = beta and NU = mu
C                           CHAR : MU = a_1
C                           SHMA : not used
C       N_MAX      :  number of moments computed
C       ORTH_POL   :  type of orthogonal polynomials used
C
C                   ='CHEB' : Chebychev polynomials
C                   ='KRAW' : Krawtchouk polynomials
C                   ='HAHN' : Hahn polynomials
C                   ='MEIX' : Meixner polynomials
C                   ='CHAR' : Charlier polynomials
C                   ='SHMA' : Shmaliy polynomials
C       W          :  w coefficient
C       DN2        :  dn square coefficient
C
C
C  Output parameters:
C
C       ON         :  discrete orthogonal polynomial 
C
C
C   Author:  D. Sébilleau 
C
C                                          Last modified :  2 May 2014
C
C
      PARAMETER (NMAX=9999)
C
      REAL*8 ON(0:NMAX,0:N_P-1),X_P,MU,NU,XN,XN1
      REAL*8 W(0:N_P-1),DN2(0:N_MAX)
      REAL*8 FACTLN,GAMMLN,COEF,CS
C
      CHARACTER*4 ORTH_POL
C
      X_P=DBLE(N_P)
      CS=-1.D0
C
      DO N=0,N_MAX-1
C
        XN=DFLOAT(N)
        XN1=XN+1.D0
        CS=-CS
C
        IF(ORTH_POL.EQ.'CHEB') THEN
C
C  Discrete Chebyshev polynomials
C
          ON(N,0)=CS*DEXP(GAMMLN(X_P)-GAMMLN(X_P-XN))*DSQRT(W(0)/DN2(N))
          COEF=(1.D0+XN*XN1/(1.D0-X_P))
          ON(N,1)=COEF*DSQRT(W(1)/W(0))*ON(N,0)
C
        ELSEIF(ORTH_POL.EQ.'KRAW') THEN
C
C  Discrete Krawtchouk polynomials
C
           ON(N,0)=DSQRT(DEXP(N*DLOG(MU)+(X_P-N)*DLOG(1.D0-MU)+
     1                        FACTLN(N_P)-FACTLN(N)-FACTLN(N_P-N)))
           ON(N,1)=ON(N,0)*(X_P*MU-N)/DSQRT(X_P*MU*(1.D0-MU))
C
        ELSEIF(ORTH_POL.EQ.'HAHN') THEN
C
C  Discrete Hahn polynomials
C
          ON(N,0)=CS*DEXP(GAMMLN(X_P)-GAMMLN(X_P-XN)+GAMMLN(XN1+NU)-
     1                 GAMMLN(XN1)-GAMMLN(NU+1.D0))*DSQRT(W(0)/DN2(N))
          COEF=((XN1+NU)*(X_P-XN1)-XN*(X_P+MU-1.D0))/
     1         ((NU+1.D0)*(X_P-1.D0))
          ON(N,1)=COEF*DSQRT(W(1)/W(0))*ON(N,0)
C
        ELSEIF(ORTH_POL.EQ.'MEIX') THEN
C
C  Discrete Meixner polynomials
C
          ON(N,0)=DEXP(GAMMLN(MU+XN)-GAMMLN(MU))*DSQRT(W(0)/DN2(N))
          COEF=(NU*(XN+MU)-XN)/(NU*MU)
          ON(N,1)=COEF*DSQRT(W(1)/W(0))*ON(N,0)
C
        ELSEIF(ORTH_POL.EQ.'CHAR') THEN
C
C  Discrete Charlier polynomials
C
          ON(N,0)=DSQRT(W(0)/DN2(N))
          COEF=(MU-XN)/MU
          ON(N,1)=COEF*DSQRT(W(1)/W(0))*ON(N,0)
C
        ENDIF
      ENDDO
C
      END
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION FACTLN(N)  
C
C            THIS SUBROUTINE IS TAKEN FROM THE BOOK :
C          "NUMERICAL RECIPES : THE ART OF SCIENTIFIC
C           COMPUTING" PAR W.H. PRESS, B.P. FLANNERY,
C               S.A. TEUKOLSKY ET W.T. VETTERLING 
C               (CAMBRIDGE UNIVERSITY PRESS 1992)
C
C                       p. 208
C
C     Logarithm of factorial function
C
      INTEGER N 
C
      REAL*8 A(100),GAMMLN  
C
      SAVE A  
C
      DATA A/100*-1.D0/  
C
      IF (N.LT.0) THEN
        PRINT *, 'NEGATIVE FACTORIAL IN FACTLN : N = ', N  
        STOP
      END IF
C
      IF (N.LE.99) THEN  
        IF (A(N+1).LT.0.D0) A(N+1)=GAMMLN(N+1.D0)  
        FACTLN=A(N+1)  
      ELSE  
        FACTLN=GAMMLN(N+1.D0)  
      ENDIF  
C
      RETURN 
C
      END  
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION GAMMLN(XX)  
C
C            THIS SUBROUTINE IS TAKEN FROM THE BOOK :
C          "NUMERICAL RECIPES : THE ART OF SCIENTIFIC
C           COMPUTING" PAR W.H. PRESS, B.P. FLANNERY,
C               S.A. TEUKOLSKY ET W.T. VETTERLING 
C               (CAMBRIDGE UNIVERSITY PRESS 1992)
C
C                       p. 207
C
C     Logarithm of Gamma function
C
      INTEGER J  
C
      REAL*8 XX 
C 
      DOUBLE PRECISION SER,STP,TMP,X,Y,COF(6)
C  
      SAVE COF,STP  
C
      DATA COF,STP/76.18009172947146D0,-86.50532032941677D0,  
     1             24.01409824083091D0,-1.231739572450155D0,  
     2               .1208650973866179D-2,-.5395239384953D-5,
     3                2.5066282746310005D0/ 
C 
      X=XX  
      Y=X  
      TMP=X+5.5D0  
      TMP=(X+0.5D0)*LOG(TMP)-TMP  
      SER=1.000000000190015D0  
C
      DO 11 J=1,6  
        Y=Y+1.D0  
        SER=SER+COF(J)/Y  
11    CONTINUE 
C 
      GAMMLN=TMP+DLOG(STP*SER/X)  
C
      RETURN  
C
      END  
C
C=======================================================================
C
      SUBROUTINE GAULEG(X1,X2,X,W,N)
C
C            This subroutine is taken from the book :
C          "Numerical Recipes : The Art of Scientific
C           Computing" par W.H. PRESS, B.P. FLANNERY,
C               S.A. TEUKOLSKY et W.T. VETTERLING 
C               (Cambridge University Press 1992)
C
C                       p. 145
C
C
C  Input parameters:
C
C       X1       :  lower limit of integration
C       X2       :  upper limit of integration 
C       N        :  order of the Gauss-Legendre quadrature formula
C
C
C  Output parameters:
C
C       X        :  abscissas for Gauss-Legendre N-point quadrature formula
C       W        :  weights for Gauss-Legendre N-point quadrature formula
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*4 X1,X2,X(N),W(N)
C
      PARAMETER (EPS=3.D-14)
C
      M=(N+1)/2 
      XM=0.5D0*DBLE(X2+X1) 
      XL=0.5D0*DBLE(X2-X1) 
C
      DO 12 I=1,M 
C
        Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0)) 
1       CONTINUE 
          P1=1.D0 
          P2=0.D0 
          DO 11 J=1,N 
            P3=P2 
            P2=P1 
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J 
11        CONTINUE 
          PP=N*(Z*P1-P2)/(Z*Z-1.D0) 
          Z1=Z 
          Z=Z1-P1/PP 
        IF(ABS(Z-Z1).GT.EPS)GO TO 1 
        X(I)=REAL(XM-XL*Z)
        X(N+1-I)=REAL(XM+XL*Z)
        W(I)=REAL(2.D0*XL/((1.D0-Z*Z)*PP*PP))
        W(N+1-I)=W(I) 
C
12    CONTINUE 
C
      RETURN 
C
      END
C
C=======================================================================
C
      SUBROUTINE POLLEG(NC,X,PL)                                        
C
C  This routine computes the Legendre polynomials up to order NC-1
C
C
C  Input parameters:
C
C       NC       :  number of l values to be computed (from 0 to NC-1)
C       X        :  argument of the Legendre polynomial
C
C
C  Output parameters:
C
C       PL       :  Legendre polynomials
C
C
      PARAMETER(NMAX=9999)
C
      DIMENSION PL(0:NMAX) 
C
      PL(0)=1.                                                          
      PL(1)=X                                                           
      DO 10 L=2,NC-1                                                    
        L1=L-1
        L2=L-2                                                            
        L3=2*L-1                                                          
        PL(L)=(X*FLOAT(L3)*PL(L1)-FLOAT(L1)*PL(L2))/FLOAT(L)
  10  CONTINUE
C
      RETURN                                                            
C
      END                                                               
C
C=======================================================================
C
      FUNCTION HFBICO(N,K)
C
C  Binomial coefficient for half integer
C
C   Author :  K. Hatada  
C
C                                          Last modified : 21 Jan 2013
      INTEGER K,N
C
      REAL*4 HFBICO
      REAL*4 FACTLNS
C
      DATA PI /3.141593/
C
      IF (0.5*REAL(N)+1.0.GE.REAL(K)) THEN
        HFBICO=(EXP(GAMMLNS(REAL(N)*0.5+1.0)-FACTLNS(K)-
     1             GAMMLNS(REAL(N)*0.5-REAL(K)+1.0)))
      ELSE
        HFBICO=(EXP(GAMMLNS(REAL(N)*0.5+1.0)-FACTLNS(K)+
     1             GAMMLNS(REAL(K)-REAL(N)*0.5-1.0)))*
     2         (REAL(K)-REAL(N)*0.5-1.0)*
     3         SIN((REAL(N)*0.5-REAL(K)+1.0)*PI)/PI
      END IF
C     
      END
C
C=======================================================================
C
      FUNCTION BICO(N,K)
C
C            THIS SUBROUTINE IS TAKEN FROM THE BOOK :
C          "NUMERICAL RECIPES : THE ART OF SCIENTIFIC
C           COMPUTING" PAR W.H. PRESS, B.P. FLANNERY,
C               S.A. TEUKOLSKY ET W.T. VETTERLING 
C               (CAMBRIDGE UNIVERSITY PRESS 1992)
C
C                       p. 208
C
C     BINOMIAL COEFFICIENT
C
      INTEGER K,N
C
      REAL*4 BICO
C
      REAL*4 FACTLNS
C
      BICO=NINT(EXP(FACTLNS(N)-FACTLNS(K)-FACTLNS(N-K)))
C
      RETURN
C
C    THE NEAREST-INTEGER FUNCTION CLEANS UP ROUNDOFF ERROR 
C             FOR SMALLER VALUES OF N AND K.
C
      END
C
C=======================================================================
C
      FUNCTION FACTLNS(N)  
C
C            THIS SUBROUTINE IS TAKEN FROM THE BOOK :
C          "NUMERICAL RECIPES : THE ART OF SCIENTIFIC
C           COMPUTING" PAR W.H. PRESS, B.P. FLANNERY,
C               S.A. TEUKOLSKY ET W.T. VETTERLING 
C               (CAMBRIDGE UNIVERSITY PRESS 1992)
C
C                       p. 208
C
C     FACTORIAL (single precision)
C
      INTEGER N 
C
      REAL*4 FACTLNS
      REAL*4 A(100),GAMMLNS 
C
      SAVE A  
C
      DATA A/100*-1./  
C
      IF (N.LT.0) THEN
        PRINT *, 'NEGATIVE FACTORIAL IN FACTLNS : N = ', N  
        STOP
      END IF
C
      IF (N.LE.99) THEN  
        IF (A(N+1).LT.0.) A(N+1)=GAMMLNS(N+1.)  
        FACTLNS=A(N+1)  
      ELSE  
        FACTLNS=GAMMLNS(N+1.)  
      ENDIF  
C
      RETURN 
C
      END  
C
C=======================================================================
C
      FUNCTION GAMMLNS(XX)  
C
C            THIS SUBROUTINE IS TAKEN FROM THE BOOK :
C          "NUMERICAL RECIPES : THE ART OF SCIENTIFIC
C           COMPUTING" PAR W.H. PRESS, B.P. FLANNERY,
C               S.A. TEUKOLSKY ET W.T. VETTERLING 
C               (CAMBRIDGE UNIVERSITY PRESS 1992)
C
C                       p. 207
C
C     Logarithm of Gamma function (single precision)
C
      INTEGER J  
C
      REAL*4 GAMMLNS,XX 
C 
      DOUBLE PRECISION SER,STP,TMP,X,Y,COF(6)
C  
      SAVE COF,STP  
C
      DATA COF,STP/76.18009172947146D0,-86.50532032941677D0,  
     1             24.01409824083091D0,-1.231739572450155D0,  
     2               .1208650973866179D-2,-.5395239384953D-5,
     3                2.5066282746310005D0/ 
C 
      X=DBLE(XX  )
      Y=X  
      TMP=X+5.5D0  
      TMP=(X+0.5D0)*LOG(TMP)-TMP  
      SER=1.000000000190015D0  
C
      DO 11 J=1,6  
        Y=Y+1.D0  
        SER=SER+COF(J)/Y  
11    CONTINUE 
C 
      GAMMLNS=REAL(TMP+LOG(STP*SER/X))
C
      RETURN  
C
      END  
C
C=======================================================================
C
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
C
C  This subroutine constructs a cubic spline 
C
C    Taken from "Numerical Recipes in Fortran 77 second edition"
C
C      from W. H. Press, S. A. Teukolsky, W. T. Vetterling 
C      and B. P. Flannery, p. 109
C
C
C  Input parameters:
C
C       X        :  x coordinates of the input function
C       Y        :  y coordinates of the input function  
C       N        :  number of points of the input function
C       YP1      :  value of first derivative of interpolating function at point 1
C       YPN      :  value of first derivative of interpolating function at point N
C
C
C  Output parameters:
C
C       Y2       :  y coordinates second derivative of the interpolating function
C
C
      PARAMETER(NMAX=1000)
C
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
C
      IF(YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
C
      DO I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     1       /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      ENDDO
C
      IF(YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
C
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
C
      DO K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
      ENDDO
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y,*)
C
C  This subroutine performs a cubic spline interpolation
C
C    Taken from "Numerical Recipes in Fortran 77 second edition"
C
C      from W. H. Press, S. A. Teukolsky, W. T. Vetterling 
C      and B. P. Flannery, p. 110
C
C
C  Input parameters:
C
C       XA       :  x coordinates of the input function
C       YA       :  y coordinates of the input function  
C       Y2A      :  y coordinates second derivative of the interpolating function
C                      (output of subroutine SPLINE)
C       N        :  number of points of the input function
C       X        :  x value at which interpolation is made
C
C
C  Output parameters:
C
C       Y        :  cubic-spline interpolated value
C
C
      DIMENSION XA(N),YA(N),Y2A(N)
C
      KLO=1
      KHI=N
C
   1  IF(KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X) THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
        GOTO 1
      ENDIF
C
      H=XA(KHI)-XA(KLO)
C
      IF(H.EQ.0.) RETURN 1
C
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     1         ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE LOCATE(XX,N,X,J,K)
C
C
C            This subroutine is taken from the book :
C          "Numerical Recipes : The Art of Scientific
C           Computing" par W.H. PRESS, B.P. FLANNERY,
C               S.A. TEUKOLSKY et W.T. VETTERLING 
C               (Cambridge University Press 1992)
C
C    It performs a search in an ordered table using a bisection method.
C    Given a monotonic array XX(1:N) and a value X, it returns J such 
C                  that X is between XX(J) and XX(J+1). 
C
C
C  Input parameters:
C
C       XX       :  coordinates of the input array
C       N        :  number of points of the input array
C       X        :  x value for which the nearest array index is sought
C
C
C  Output parameters:
C
C       J        :  index for which XX(J) <= > <= XX(J+1)
C
C
      INTEGER J,N
      INTEGER JL,JM,JU
C
      REAL X,XX(N)
C
      JL=0
      JU=N+1
  10  IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GOTO 10
      ENDIF
      IF(K.EQ.1) THEN
        J=JL+1
      ELSE
        XU=ABS(X-XX(JU))
        IF(XU.LT.0.0001) THEN
          J=JU
        ELSE
          J=JL
        ENDIF
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE DERIV(F,N,F1,F2,F3,F4,F5,N_POINTS,H,I_FLAG)
C
C   This subroutine computes the first (F1), second (F2) 
C     third (F3), fourth (F4) and fifth (F5) derivatives of function F.
C
C   The general formula used is generally a central difference formula,
C     except for the first two points (forward difference formula) 
C      and the last two or three points (backward difference formula).
C
C   N_POINTS : number of points used for the calculation
C                 (2 <= N_POINTS <= 6)
C
C   I_FLAG   : = 1 : only first derivative computed
C              = 2 : first and second derivative computed
C              = 3 : idem 2 + third derivative computed
C              = 4 : idem 3 + fourth derivative computed
C              = 5 : idem 4 + fifth derivative computed
C
C   H        : step
C
C   References : A. K. Singh and G. R. Thorpe, 
C                RGMIA Res. Rep. Coll., 2(6), Article 7, 1999.
C
C                T. F. Guidry, 
C               http://www.trentfguidry.net/post/2010/09/04/Numerical-differentiation-formulas.aspx
C
C
C  Input parameters:
C
C       F        :  y coordinates of the input file
C       N        :  dimension of the arrays
C       N_POINTS :  number of points in the file
C       H        :  step of the input file
C       I_FLAG   :  flag to select the number of derivatives computed
C
C
C  Output parameters:
C
C       Fn       :  order n derivative of F
C
C
C   Author :  D. Sébilleau
C                                          Last version : 22 Dec 2014
C 
C
      REAL*4 F(N),F1(N),F2(N),F3(N),F4(N)
      REAL*4 F5(N)
      REAL*4 A(10,0:10),B(10,-10:0),C(10,-10:10)
      REAL*4 STEP1,STEP2,STEP3,STEP4,STEP5
C
      INTEGER N_MIN(10)
C
      DATA N_MIN /1,2,3,6,6,10,0,0,0,0/
C
C   Check of the consistency of the number of points
C         and the order of the derivative
C
      IF(N_POINTS.LT.(I_FLAG+1)) THEN
        N_POINTS=I_FLAG+1
      ENDIF
C
C  Check if the number of points N in function F is sufficient 
C    to computes the derivatives in view of the algorithm used
C
  20  IF(N.LT.N_MIN(N_POINTS)) THEN
        N_POINTS=N_POINTS-1
        WRITE(6,10) N_POINTS
        GOTO 20
      ENDIF
C
C  Computation of the derivative(s) using a N-POINTS formula
C
      IF(N_POINTS.EQ.2) THEN
C
C            2-POINT FORMULA :
C
C            no second derivative
C
        STEP1=H
C
        CALL COEF_DERIV(N_POINTS,A,B,C)
C
        F1(1)=(A(1,0)*F(1)+A(1,1)*F(2))/STEP1
C
        DO JP=2,N
          F1(JP)=(B(1,0)*F(JP)+B(1,-1)*F(JP-1))/STEP1
        ENDDO
C
      ELSEIF(N_POINTS.EQ.3) THEN
C
C            3-POINT FORMULA :
C
C            no third derivative
C
        STEP1=2.*H
        STEP2=H*H
C
        CALL COEF_DERIV(N_POINTS,A,B,C)
C
        F1(1)=(A(1,0)*F(1)+A(1,1)*F(2)+A(1,2)*F(3))/STEP1
C
        F1(N)=(B(1,0)*F(N)+B(1,-1)*F(N-1)+B(1,-2)*F(N-2))/STEP1
C
        IF(I_FLAG.GE.2) THEN
          F2(1)=(A(2,0)*F(1)+A(2,1)*F(2)+A(2,2)*F(3))/STEP2
C
          F2(N)=(B(2,0)*F(N)+B(2,-1)*F(N-1)+B(2,-2)*F(N-2))/STEP2
        ENDIF
C
        DO JP=2,N-1 
          F1(JP)=(C(1,-1)*F(JP-1)+C(1,0)*F(JP)+C(1,1)*F(JP+1))/STEP1
C
          IF(I_FLAG.GE.2) THEN
            F2(JP)=(C(2,-1)*F(JP-1)+C(2,0)*F(JP)+C(2,1)*F(JP+1))/STEP2
          ENDIF
C
        ENDDO       
C
      ELSEIF(N_POINTS.EQ.4) THEN
C
C            4-POINT FORMULA :
C
C            no fourth derivative
C
        STEP1=6.*H
        STEP2=H*H
        STEP3=STEP2*H
C
        CALL COEF_DERIV(N_POINTS,A,B,C)
C
        F1(1)=(A(1,0)*F(1)+A(1,1)*F(2)+A(1,2)*F(3)+A(1,3)*F(4))/STEP1
        F1(2)=(A(1,0)*F(2)+A(1,1)*F(3)+A(1,2)*F(4)+A(1,3)*F(5))/STEP1
        F1(3)=(A(1,0)*F(3)+A(1,1)*F(4)+A(1,2)*F(5)+A(1,3)*F(6))/STEP1
C
        IF(I_FLAG.GE.2) THEN
          F2(1)=(A(2,0)*F(1)+A(2,1)*F(2)+A(2,2)*F(3)+A(2,3)*F(4))/STEP2
          F2(2)=(A(2,0)*F(2)+A(2,1)*F(3)+A(2,2)*F(4)+A(2,3)*F(5))/STEP2
          F2(3)=(A(2,0)*F(3)+A(2,1)*F(4)+A(2,2)*F(5)+A(2,3)*F(6))/STEP2
        ENDIF
C
        IF(I_FLAG.GE.3) THEN
          F3(1)=(A(3,0)*F(1)+A(3,1)*F(2)+A(3,2)*F(3)+A(3,3)*F(4))/STEP3
          F3(2)=(A(3,0)*F(2)+A(3,1)*F(3)+A(3,2)*F(4)+A(3,3)*F(5))/STEP3
          F3(3)=(A(3,0)*F(3)+A(3,1)*F(4)+A(3,2)*F(5)+A(3,3)*F(6))/STEP3
        ENDIF
C
        DO JP=4,N
          F1(JP)=(B(1,0)*F(JP)+B(1,-1)*F(JP-1)+B(1,-2)*F(JP-2)+
     1            B(1,-3)*F(JP-3))/STEP1
C
          IF(I_FLAG.GE.2) THEN
            F2(JP)=(B(2,0)*F(JP)+B(2,-1)*F(JP-1)+B(2,-2)*F(JP-2)+
     1              B(2,-3)*F(JP-3))/STEP2
          ENDIF
C
          IF(I_FLAG.GE.3) THEN
            F3(JP)=(B(3,0)*F(JP)+B(3,-1)*F(JP-1)+B(3,-2)*F(JP-2)+
     1              B(3,-3)*F(JP-3))/STEP3
          ENDIF
C
        ENDDO
C
      ELSEIF(N_POINTS.EQ.5) THEN
C
C            5-POINT FORMULA :
C
C            no fifth derivative
C
        STEP1=12.*H
        STEP2=12.*H*H
        STEP3=2.*H*H*H
        STEP4=H*H*H*H
C
        CALL COEF_DERIV(N_POINTS,A,B,C)
C
        F1(1)=(A(1,0)*F(1)+A(1,1)*F(2)+A(1,2)*F(3)+A(1,3)*F(4)+
     1         A(1,4)*F(5))/STEP1
        F1(2)=(A(1,0)*F(2)+A(1,1)*F(3)+A(1,2)*F(4)+A(1,3)*F(5)+
     1         A(1,4)*F(6))/STEP1
C
        F1(N-1)=(B(1,0)*F(N-1)+B(1,-1)*F(N-2)+B(1,-2)*F(N-3)+
     1           B(1,-3)*F(N-4)+B(1,-4)*F(N-5))/STEP1
        F1(N)=(B(1,0)*F(N)+B(1,-1)*F(N-1)+B(1,-2)*F(N-2)+
     1         B(1,-3)*F(N-3)+B(1,-4)*F(N-4))/STEP1
C
        IF(I_FLAG.GE.2) THEN
          F2(1)=(A(2,0)*F(1)+A(2,1)*F(2)+A(2,2)*F(3)+A(2,3)*F(4)+
     1           A(2,4)*F(5))/STEP2
          F2(2)=(A(2,0)*F(2)+A(2,1)*F(3)+A(2,2)*F(4)+A(2,3)*F(5)+
     1           A(2,4)*F(6))/STEP2
C
          F2(N-1)=(B(2,0)*F(N-1)+B(2,-1)*F(N-2)+B(2,-2)*F(N-3)+
     1             B(2,-3)*F(N-4)+B(2,-4)*F(N-5))/STEP2
          F2(N)=(B(2,0)*F(N)+B(2,-1)*F(N-1)+B(2,-2)*F(N-2)+
     1           B(2,-3)*F(N-3)+B(2,-4)*F(N-4))/STEP2
        ENDIF
C
        IF(I_FLAG.GE.3) THEN
          F3(1)=(A(3,0)*F(1)+A(3,1)*F(2)+A(3,2)*F(3)+A(3,3)*F(4)+
     1           A(3,4)*F(5))/STEP3
          F3(2)=(A(3,0)*F(2)+A(3,1)*F(3)+A(3,2)*F(4)+A(3,3)*F(5)+
     1           A(3,4)*F(6))/STEP3
C
          F3(N-1)=(B(3,0)*F(N-1)+B(3,-1)*F(N-2)+B(3,-2)*F(N-3)+
     1             B(3,-3)*F(N-4)+B(3,-4)*F(N-5))/STEP3
          F3(N)=(B(3,0)*F(N)+B(3,-1)*F(N-1)+B(3,-2)*F(N-2)+
     1           B(3,-3)*F(N-3)+B(3,-4)*F(N-4))/STEP3
        ENDIF
C
        IF(I_FLAG.GE.4) THEN
          F4(1)=(A(4,0)*F(1)+A(4,1)*F(2)+A(4,2)*F(3)+A(4,3)*F(4)+
     1           A(4,4)*F(5))/STEP4
          F4(2)=(A(4,0)*F(2)+A(4,1)*F(3)+A(4,2)*F(4)+A(4,3)*F(5)+
     1           A(4,4)*F(6))/STEP4
C
          F4(N-1)=(B(4,0)*F(N-1)+B(4,-1)*F(N-2)+B(4,-2)*F(N-3)+
     1             B(4,-3)*F(N-4)+B(4,-4)*F(N-5))/STEP4
          F4(N)=(B(4,0)*F(N)+B(4,-1)*F(N-1)+B(4,-2)*F(N-2)+
     1           B(4,-3)*F(N-3)+B(4,-4)*F(N-4))/STEP4
        ENDIF
C
        DO JP=3,N-2
C
          F1(JP)=(C(1,-2)*F(JP-2)+C(1,-1)*F(JP-1)+C(1,0)*F(JP)+
     1            C(1,1)*F(JP+1)+C(1,2)*F(JP+2))/STEP1
C
          IF(I_FLAG.GE.2) THEN
            F2(JP)=(C(2,-2)*F(JP-2)+C(2,-1)*F(JP-1)+C(2,0)*F(JP)+
     1              C(2,1)*F(JP+1)+C(2,2)*F(JP+2))/STEP2
          ENDIF
C
          IF(I_FLAG.GE.3) THEN
            F3(JP)=(C(3,-2)*F(JP-2)+C(3,-1)*F(JP-1)+C(3,0)*F(JP)+
     1              C(3,1)*F(JP+1)+C(3,2)*F(JP+2))/STEP3
          ENDIF
C
          IF(I_FLAG.GE.4) THEN
            F4(JP)=(C(4,-2)*F(JP-2)+C(4,-1)*F(JP-1)+C(4,0)*F(JP)+
     1              C(4,1)*F(JP+1)+C(4,2)*F(JP+2))/STEP4
          ENDIF
C
        ENDDO
C
      ELSEIF(N_POINTS.EQ.6) THEN
C
C            6-POINT FORMULA :
C
        STEP1=60.*H
        STEP2=12.*H*H
        STEP3=4.*H*H*H
        STEP4=H*H*H*H
        STEP5=STEP4*H
C
        CALL COEF_DERIV(N_POINTS,A,B,C)
C
        F1(1)=(A(1,0)*F(1)+A(1,1)*F(2)+A(1,2)*F(3)+A(1,3)*F(4)+
     1         A(1,4)*F(5)+A(1,5)*F(6))/STEP1
        F1(2)=(A(1,0)*F(2)+A(1,1)*F(3)+A(1,2)*F(4)+A(1,3)*F(5)+
     1         A(1,4)*F(6)+A(1,5)*F(7))/STEP1
        F1(3)=(A(1,0)*F(3)+A(1,1)*F(4)+A(1,2)*F(5)+A(1,3)*F(6)+
     1         A(1,4)*F(7)+A(1,5)*F(8))/STEP1
        F1(4)=(A(1,0)*F(4)+A(1,1)*F(5)+A(1,2)*F(6)+A(1,3)*F(7)+
     1         A(1,4)*F(8)+A(1,5)*F(9))/STEP1
        F1(5)=(A(1,0)*F(5)+A(1,1)*F(6)+A(1,2)*F(7)+A(1,3)*F(8)+
     1         A(1,4)*F(9)+A(1,5)*F(10))/STEP1
C
        IF(I_FLAG.GE.2) THEN
          F2(1)=(A(2,0)*F(1)+A(2,1)*F(2)+A(2,2)*F(3)+A(2,3)*F(4)+
     1           A(2,4)*F(5)+A(2,5)*F(6))/STEP2
          F2(2)=(A(2,0)*F(2)+A(2,1)*F(3)+A(2,2)*F(4)+A(2,3)*F(5)+
     1           A(2,4)*F(6)+A(2,5)*F(7))/STEP2
          F2(3)=(A(2,0)*F(3)+A(2,1)*F(4)+A(2,2)*F(5)+A(2,3)*F(6)+
     1           A(2,4)*F(7)+A(2,5)*F(8))/STEP2
          F2(4)=(A(2,0)*F(4)+A(2,1)*F(5)+A(2,2)*F(6)+A(2,3)*F(7)+
     1           A(2,4)*F(8)+A(2,5)*F(9))/STEP2
          F2(5)=(A(2,0)*F(5)+A(2,1)*F(6)+A(2,2)*F(7)+A(2,3)*F(8)+
     1           A(2,4)*F(9)+A(2,5)*F(10))/STEP2
        ENDIF
C
        IF(I_FLAG.GE.3) THEN
          F3(1)=(A(3,0)*F(1)+A(3,1)*F(2)+A(3,2)*F(3)+A(3,3)*F(4)+
     1           A(3,4)*F(5)+A(3,5)*F(6))/STEP3
          F3(2)=(A(3,0)*F(2)+A(3,1)*F(3)+A(3,2)*F(4)+A(3,3)*F(5)+
     1           A(3,4)*F(6)+A(3,5)*F(7))/STEP3
          F3(3)=(A(3,0)*F(3)+A(3,1)*F(4)+A(3,2)*F(5)+A(3,3)*F(6)+
     1           A(3,4)*F(7)+A(3,5)*F(8))/STEP3
          F3(4)=(A(3,0)*F(4)+A(3,1)*F(5)+A(3,2)*F(6)+A(3,3)*F(7)+
     1           A(3,4)*F(8)+A(3,5)*F(9))/STEP3
          F3(5)=(A(3,0)*F(5)+A(3,1)*F(6)+A(3,2)*F(7)+A(3,3)*F(8)+
     1           A(3,4)*F(9)+A(3,5)*F(10))/STEP3
        ENDIF
C
        IF(I_FLAG.GE.4) THEN
          F4(1)=(A(4,0)*F(1)+A(4,1)*F(2)+A(4,2)*F(3)+A(4,3)*F(4)+
     1           A(4,4)*F(5)+A(4,5)*F(6))/STEP4
          F4(2)=(A(4,0)*F(2)+A(4,1)*F(3)+A(4,2)*F(4)+A(4,3)*F(5)+
     1           A(4,4)*F(6)+A(4,5)*F(7))/STEP4
          F4(3)=(A(4,0)*F(3)+A(4,1)*F(4)+A(4,2)*F(5)+A(4,3)*F(6)+
     1           A(4,4)*F(7)+A(4,5)*F(8))/STEP4
          F4(4)=(A(4,0)*F(4)+A(4,1)*F(5)+A(4,2)*F(6)+A(4,3)*F(7)+
     1           A(4,4)*F(8)+A(4,5)*F(9))/STEP4
          F4(5)=(A(4,0)*F(5)+A(4,1)*F(6)+A(4,2)*F(7)+A(4,3)*F(8)+
     1           A(4,4)*F(9)+A(4,5)*F(10))/STEP4
        ENDIF
C
        IF(I_FLAG.GE.5) THEN
          F5(1)=(A(5,0)*F(1)+A(5,1)*F(2)+A(5,2)*F(3)+A(5,3)*F(4)+
     1           A(5,4)*F(5)+A(5,5)*F(6))/STEP5
          F5(2)=(A(5,0)*F(2)+A(5,1)*F(3)+A(5,2)*F(4)+A(5,3)*F(5)+
     1           A(5,4)*F(6)+A(5,5)*F(7))/STEP5
          F5(3)=(A(5,0)*F(3)+A(5,1)*F(4)+A(5,2)*F(5)+A(5,3)*F(6)+
     1           A(5,4)*F(7)+A(5,5)*F(8))/STEP5
          F5(4)=(A(5,0)*F(4)+A(5,1)*F(5)+A(5,2)*F(6)+A(5,3)*F(7)+
     1           A(5,4)*F(8)+A(5,5)*F(9))/STEP5
          F5(5)=(A(5,0)*F(5)+A(5,1)*F(6)+A(5,2)*F(7)+A(5,3)*F(8)+
     1           A(5,4)*F(9)+A(5,5)*F(10))/STEP5
        ENDIF
C
        DO JP=6,N
C
          F1(JP)=(B(1,0)*F(JP)+B(1,-1)*F(JP-1)+B(1,-2)*F(JP-2)+
     1            B(1,-3)*F(JP-3)+B(1,-4)*F(JP-4)+
     2            B(1,-5)*F(JP-5))/STEP1
C
          IF(I_FLAG.GE.2) THEN
            F2(JP)=(B(2,0)*F(JP)+B(2,-1)*F(JP-1)+B(2,-2)*F(JP-2)+
     1              B(2,-3)*F(JP-3)+B(2,-4)*F(JP-4)+
     2              B(2,-5)*F(JP-5))/STEP2
          ENDIF
C
          IF(I_FLAG.GE.3) THEN
            F3(JP)=(B(3,0)*F(JP)+B(3,-1)*F(JP-1)+B(3,-2)*F(JP-2)+
     1              B(3,-3)*F(JP-3)+B(3,-4)*F(JP-4)+
     2              B(3,-5)*F(JP-5))/STEP3
          ENDIF
C
          IF(I_FLAG.GE.4) THEN
            F4(JP)=(B(4,0)*F(JP)+B(4,-1)*F(JP-1)+B(4,-2)*F(JP-2)+
     1              B(4,-3)*F(JP-3)+B(4,-4)*F(JP-4)+
     2              B(4,-5)*F(JP-5))/STEP4
          ENDIF
C
          IF(I_FLAG.GE.5) THEN
            F5(JP)=(B(5,0)*F(JP)+B(5,-1)*F(JP-1)+B(5,-2)*F(JP-2)+
     1              B(5,-3)*F(JP-3)+B(5,-4)*F(JP-4)+
     2              B(5,-5)*F(JP-5))/STEP5
          ENDIF
C
        ENDDO
C
      ENDIF
C
C  Format
C
  10  FORMAT(//,10X,'<<<<<  NOT ENOUGH POINTS IN FUNCTION:  >>>>>',
     1        /,10X,'<<<<<   USING ',I1,'-POINTS FORMULA    >>>>>',//)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE COEF_DERIV(NP,A,B,C)
C
C      This subroutine computes the coefficients for the 
C             NP-point derivation with 1 < NP < 8 
C
C      Derivatives up to order (NP-1) can be computed from 
C            these coefficients (limited to order 5)
C
C  Input parameters:
C
C       * NP       : number of points of the derivation

C
C  Output parameters:
C
C       * A(ND,NP) : coefficients of the derivation for the forward 
C                        difference scheme
C       * B(ND,NP) : coefficients of the derivation for the backward 
C                        difference scheme
C       * C(ND,NP) : coefficients of the derivation for the central
C                        difference scheme
C
C               with ND the order of the derivation
C
C  References:  T. F. Guidry, 
C               http://www.trentfguidry.net/post/2010/09/04/Numerical-differentiation-formulas.aspx


C  Note: the coefficients are computed for three different schemes:
C
C               = F : forward difference 
C               = B : backward difference 
C               = C : central difference (Stirling)
C
C        The order of the coefficients is the following:
C
C               = F :  A(0)*F(I) + A(1)*F(I+1) + ...
C               = B :  B(0)*F(I) + B(-1)*F(I-1) + ...
C               = C :  ... + C(-1)*F(I-1) + C(0)*F(I) + C(1)*F(I+1) + ...
C
C                                                         
C   Author :  D. Sébilleau
C
C                                          Last modified : 19 Nov 2014
C   
      INTEGER NP
C
      REAL*4 A(10,0:10),B(10,-10:0),C(10,-10:10)
C
C  Initializations
C
      DO J=1,10
        DO K=0,10
          A(J,K)=0.0
        ENDDO
        DO K=-10,0
          B(J,K)=0.0
        ENDDO
        DO K=-10,10
          C(J,K)=0.0
        ENDDO
      ENDDO
C
      IF(NP.EQ.2) THEN
C
C        Forward difference scheme
C
          A(1,0)=-1.0
          A(1,1)=1.0
C
C        Backward difference scheme
C
          B(1,0)=1.0
          B(1,-1)=-1.0
C
      ELSEIF(NP.EQ.3) THEN
C
C        Forward difference scheme
C
          A(1,0)=-3.0
          A(1,1)=4.0
          A(1,2)=-1.0
C
          A(2,0)=1.0
          A(2,1)=-2.0
          A(2,2)=1.0
C
C        Backward difference scheme
C
          B(1,0)=3.0
          B(1,-1)=-4.0
          B(1,-2)=1.0
C
          B(2,0)=1.0
          B(2,-1)=-2.0
          B(2,-2)=1.0
C
C        Central difference scheme
C
          C(1,-1)=-1.0
          C(1,0)=0.0
          C(1,1)=1.0
C
          C(2,-1)=1.0
          C(2,0)=-2.0
          C(2,1)=1.0
C
      ELSEIF(NP.EQ.4) THEN
C
C        Forward difference scheme
C
          A(1,0)=-11.0
          A(1,1)=18.0
          A(1,2)=-9.0
          A(1,3)=2.0
C
          A(2,0)=2.0
          A(2,1)=-5.0
          A(2,2)=4.0
          A(2,3)=-1.0
C
          A(3,0)=-1.0
          A(3,1)=3.0
          A(3,2)=-3.0
          A(3,3)=1.0
C
C        Backward difference scheme
C
          B(1,0)=11.0
          B(1,-1)=-18.0
          B(1,-2)=9.0
          B(1,-3)=-2.0
C
          B(2,0)=2.0
          B(2,-1)=-5.0
          B(2,-2)=4.0
          B(2,-3)=-1.0
C
          B(3,0)=1.0
          B(3,-1)=-3.0
          B(3,-2)=3.0
          B(3,-3)=-1.0
C
      ELSEIF(NP.EQ.5) THEN
C
C        Forward difference scheme
C
          A(1,0)=-25.0
          A(1,1)=48.0
          A(1,2)=-36.0
          A(1,3)=16.0
          A(1,4)=-3.0
C
          A(2,0)=35.0
          A(2,1)=-104.0
          A(2,2)=114.0
          A(2,3)=-56.0
          A(2,4)=11.0
C
          A(3,0)=-5.0
          A(3,1)=18.0
          A(3,2)=-24.0
          A(3,3)=14.0
          A(3,4)=-3.0
C
          A(4,0)=1.0
          A(4,1)=-4.0
          A(4,2)=6.0
          A(4,3)=-4.0
          A(4,4)=1.0
C
C        Backward difference scheme
C
          B(1,0)=25.0
          B(1,-1)=-48.0
          B(1,-2)=36.0
          B(1,-3)=-16.0
          B(1,-4)=3.0
C
          B(2,0)=35.0
          B(2,-1)=-104.0
          B(2,-2)=114.0
          B(2,-3)=-56.0
          B(2,-4)=11.0
C
          B(3,0)=5.0
          B(3,-1)=-18.0
          B(3,-2)=24.0
          B(3,-3)=-14.0
          B(3,-4)=3.0
C
          B(4,0)=1.0
          B(4,-1)=-4.0
          B(4,-2)=6.0
          B(4,-3)=-4.0
          B(4,-4)=1.0
C
C        Central difference scheme
C
          C(1,-2)=1.0
          C(1,-1)=-8.0
          C(1,0)=0.0
          C(1,1)=8.0
          C(1,2)=-1.0
C
          C(2,-2)=-1.0
          C(2,-1)=16.0
          C(2,0)=-30.0
          C(2,1)=16.0
          C(2,2)=-1.0
C
          C(3,-2)=-1.0
          C(3,-1)=2.0
          C(3,0)=0.0
          C(3,1)=-2.0
          C(3,2)=1.0
C
          C(4,-2)=1.0
          C(4,-1)=-4.0
          C(4,0)=6.0
          C(4,1)=-4.0
          C(4,2)=1.0
C
      ELSEIF(NP.EQ.6) THEN
C
C        Forward difference scheme
C
          A(1,0)=-137.0
          A(1,1)=300.0
          A(1,2)=-300.0
          A(1,3)=200.0
          A(1,4)=-75.0
          A(1,5)=12.0
C
          A(2,0)=45.0
          A(2,1)=-154.0
          A(2,2)=214.0
          A(2,3)=-156.0
          A(2,4)=61.0
          A(2,5)=-10.0
C
          A(3,0)=-17.0
          A(3,1)=71.0
          A(3,2)=-118.0
          A(3,3)=98.0
          A(3,4)=-41.0
          A(3,5)=7.0
C
          A(4,0)=3.0
          A(4,1)=-14.0
          A(4,2)=26.0
          A(4,3)=-24.0
          A(4,4)=11.0
          A(4,5)=-2.0
C
          A(5,0)=-1.0
          A(5,1)=5.0
          A(5,2)=-10.0
          A(5,3)=10.0
          A(5,4)=-5.0
          A(5,5)=1.0
C
C        Backward difference scheme
C
          B(1,0)=137.0
          B(1,-1)=-300.0
          B(1,-2)=300.0
          B(1,-3)=-200.0
          B(1,-4)=75.0
          B(1,-5)=-12.0
C
          B(2,0)=45.0
          B(2,-1)=-154.0
          B(2,-2)=214.0
          B(2,-3)=-156.0
          B(2,-4)=61.0
          B(2,-5)=-10.0
C
          B(3,0)=17.0
          B(3,-1)=-71.0
          B(3,-2)=118.0
          B(3,-3)=-98.0
          B(3,-4)=41.0
          B(3,-5)=-7.0
C
          B(4,0)=3.0
          B(4,-1)=-14.0
          B(4,-2)=26.0
          B(4,-3)=-24.0
          B(4,-4)=11.0
          B(4,-5)=-2.0
C
          B(5,0)=1.0
          B(5,-1)=-5.0
          B(5,-2)=10.0
          B(5,-3)=-10.0
          B(5,-4)=5.0
          B(5,-5)=-1.0
C
C
      ELSEIF(NP.EQ.7) THEN
C
C        Forward difference scheme
C
          A(1,0)=-147.0
          A(1,1)=360.0
          A(1,2)=-450.0
          A(1,3)=400.0
          A(1,4)=-225.0
          A(1,5)=72.0
          A(1,6)=-10.0
C
          A(2,0)=812.0
          A(2,1)=-3132.0
          A(2,2)=5265.0
          A(2,3)=-5080.0
          A(2,4)=2970.0
          A(2,5)=-972.0
          A(2,6)=137.0
C
          A(3,0)=-49.0
          A(3,1)=232.0
          A(3,2)=-461.0
          A(3,3)=496.0
          A(3,4)=-307.0
          A(3,5)=104.0
          A(3,6)=-15.0
C
          A(4,0)=35.0
          A(4,1)=-186.0
          A(4,2)=411.0
          A(4,3)=-484.0
          A(4,4)=321.0
          A(4,5)=-114.0
          A(4,6)=17.0
C
          A(5,0)=-7.0
          A(5,1)=40.0
          A(5,2)=-95.0
          A(5,3)=120.0
          A(5,4)=-85.0
          A(5,5)=32.0
          A(5,6)=-5.0
C
C        Backward difference scheme
C
          B(1,0)=147.0
          B(1,-1)=-360.0
          B(1,-2)=450.0
          B(1,-3)=-400.0
          B(1,-4)=225.0
          B(1,-5)=-72.0
          B(1,-6)=10.0
C
          B(2,0)=812.0
          B(2,-1)=-3132.0
          B(2,-2)=5265.0
          B(2,-3)=-5080.0
          B(2,-4)=2970.0
          B(2,-5)=-972.0
          B(2,-6)=137.0
C
          B(3,0)=49.0
          B(3,-1)=-232.0
          B(3,-2)=461.0
          B(3,-3)=-496.0
          B(3,-4)=307.0
          B(3,-5)=-104.0
          B(3,-6)=15.0
C
          B(4,0)=35.0
          B(4,-1)=-186.0
          B(4,-2)=411.0
          B(4,-3)=-484.0
          B(4,-4)=321.0
          B(4,-5)=-114.0
          B(4,-6)=17.0
C
          B(5,0)=7.0
          B(5,-1)=-40.0
          B(5,-2)=95.0
          B(5,-3)=-120.0
          B(5,-4)=85.0
          B(5,-5)=-32.0
          B(5,-6)=5.0
C
C        Central difference scheme
C
          C(1,-3)=-1.0
          C(1,-2)=9.0
          C(1,-1)=-45.0
          C(1,0)=0.0
          C(1,1)=45.0
          C(1,2)=-9.0
          C(1,3)=1.0
C
          C(2,-3)=2.0
          C(2,-2)=-27.0
          C(2,-1)=270.0
          C(2,0)=-490.0
          C(2,1)=270.0
          C(2,2)=-27.0
          C(2,3)=2.0
C
          C(3,-3)=1.0
          C(3,-2)=-8.0
          C(3,-1)=13.0
          C(3,0)=0.0
          C(3,1)=-13.0
          C(3,2)=8.0
          C(3,3)=-1.0
C
          C(4,-3)=-1.0
          C(4,-2)=12.0
          C(4,-1)=-39.0
          C(4,0)=56.0
          C(4,1)=-39.0
          C(4,2)=12.0
          C(4,3)=-1.0
C
          C(5,-3)=-1.0
          C(5,-2)=4.0
          C(5,-1)=-5.0
          C(5,0)=0.0
          C(5,1)=5.0
          C(5,2)=-4.0
          C(5,3)=1.0
C
      ENDIF
C
      RETURN
C
      END 
C
C=======================================================================
C
      SUBROUTINE INTEGR_I(X,F,F_1,F_3,F_5,N_BEG,N_END,N_POINTS,
     1                    METH,N_RULE,RES)
C
C      This is the driver routine that calls the subroutine that 
C          integrates a function F(X), defined over 
C          the interval [1,N_POINTS] with constant step H
C          over the interval [N_BEG,N_END]
C
C      To increase the accuracy, it computes the integral according 
C          to different schemes. There are four ways to compute the
C          integral: int[N_BEG,N_END]                                (1)
C                    int[1,N_END]-int[1,N_BEG]                       (2)
C                    int[N_BEG,N_POINTS]-int[N_END,N_POINTS]         (3)
C                    int[1,N_POINTS]-int[1,N_BEG]-int[N_END,N_POINTS](4)
C
C      Method (4) is never used as it is equivalent either to method (2) 
C          or to method (3) in terms of accuracy
C
C      This subroutine selects the method involving the larger number
C          of points, i.e. max([N_BEG,N_END],int[1,N_BEG],[N_END,N_POINTS])
C
C
C  Input parameters:
C
C       * X        : X point of function to be integrated
C       * F        : function to be integrated
C       * F_1      : first order derivative of F
C       * F_3      : third order derivative of F
C       * F_5      : fifth order derivative of F
C       * N_BEG    : starting X point for integration of F
C       * N_END    : end X point for integration of F
C       * N_POINTS : dimensioning of F (1 to N_POINTS)
C       * METH     : integration method used
C
C                     = NCQ : Newton-Cotes 
C                     = EMS : Euler-Mac Laurin summation 
C
C       * N_RULE   : number of points used in the quadrature formula
C
C                       NCQ : Newton-Cotes quadrature rule     | Accuracy
C                          
C                           --> N_RULE =  2 : trapezoidal      |  H^3
C                           --> N_RULE =  3 : Simpson 1/3      |  H^5
C                           --> N_RULE =  4 : Simpson 3/8      |  H^5
C                           --> N_RULE =  5 : Boole/Milne      |  H^7
C                           --> N_RULE =  6 : Weddle           |  H^7
C
C                       EMS : Euler-Mac Laurin summation       | Accuracy
C 
C                           --> N_RULE =  2 (uses F_1)         |  H^5
C                           --> N_RULE =  3 (uses F_1,F_3)     |  H^7      
C                           --> N_RULE =  4 (uses F_1,F_3,F_5) |  H^9      
C
C                            BN(J) is a Bernoulli number
C
C
C  Output parameters:
C
C       * RES      : result of the integration of F over the whole 
C                       interval [1,N_END]
C
C                                                         
C   Author :  D. Sébilleau
C
C                                          Last modified : 31 Oct 2014
C     
C
      REAL*4 F(N_POINTS),F_1(N_POINTS),F_3(N_POINTS),F_5(N_POINTS)
      REAL*4 X(N_POINTS)
      REAL*4 RES,RES1,RES2
C
      CHARACTER*3 METH
C
C  Checking the number of points in the integration interval  
C    with respect to that over which the function F(X) is defined
C
      N_SIZE_I=N_END-N_BEG+1
      N_SIZE_L=N_BEG
      N_SIZE_U=N_POINTS-N_END+1
C
      N_HALF=N_POINTS/2
C
      IF(N_SIZE_I.GE.N_HALF) THEN
C
C...........  Interval of integration larger than half of  ...........
C...........      the interval of definition of F(X)        ...........
C
C                       Using method (1)
C
        CALL INTEGR_INT(X,F,F_1,F_3,F_5,N_BEG,N_END,N_POINTS,
     1                  METH,N_RULE,RES)
C
      ELSE
C
C...........  Interval of integration smaller than half of  ...........
C...........      the interval of definition of F(X)        ...........
C
        IF(N_SIZE_U.GE.N_SIZE_L) THEN
C
C                       Using method (3)
C
          CALL INTEGR_INT(X,F,F_1,F_3,F_5,N_BEG,N_POINTS,N_POINTS,
     1                    METH,N_RULE,RES1)
          CALL INTEGR_INT(X,F,F_1,F_3,F_5,N_END,N_POINTS,N_POINTS,
     1                    METH,N_RULE,RES2)
C
          RES=RES1-RES2
C
        ELSE 
C
C                       Using method (2)
C
          CALL INTEGR_INT(X,F,F_1,F_3,F_5,1,N_END,N_POINTS,
     1                    METH,N_RULE,RES1)
          CALL INTEGR_INT(X,F,F_1,F_3,F_5,1,N_BEG,N_POINTS,
     1                    METH,N_RULE,RES2)
C
          RES=RES1-RES2
C
        ENDIF
C
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE INTEGR_INT(X,F,F_1,F_3,F_5,N_BEG,N_END,N_POINTS,
     1                      METH,N_RULE,RES)
C
C      This subroutine integrates the function F(X), defined over 
C          the interval [1,N_POINTS] with constant step H
C          over the interval [N_BEG,N_END].
C
C
C  Input parameters:
C
C       * X        : X point of function to be integrated
C       * F        : function to be integrated
C       * F_1      : first order derivative of F
C       * F_3      : third order derivative of F
C       * F_5      : fifth order derivative of F
C       * N_BEG    : starting X point for integration of F
C       * N_END    : end X point for integration of F
C       * N_POINTS : dimensioning of F (1 to N_POINTS)
C       * METH     : integration method used
C
C                     = NCQ : Newton-Cotes 
C                     = EMS : Euler-Mac Laurin summation 
C
C       * N_RULE   : number of points used in the quadrature formula
C
C                       NCQ : Newton-Cotes quadrature rule     | Accuracy
C                          
C                           --> N_RULE =  2 : trapezoidal      |  H^3
C                           --> N_RULE =  3 : Simpson 1/3      |  H^5
C                           --> N_RULE =  4 : Simpson 3/8      |  H^5
C                           --> N_RULE =  5 : Boole/Milne      |  H^7
C                           --> N_RULE =  6 : Weddle           |  H^7
C
C                       EMS : Euler-Mac Laurin summation       | Accuracy
C 
C                           --> N_RULE =  2 (uses F_1)         |  H^5
C                           --> N_RULE =  3 (uses F_1,F_3)     |  H^7      
C                           --> N_RULE =  4 (uses F_1,F_3,F_5) |  H^9      
C
C                            BN(J) is a Bernoulli number
C
C
C  Output parameters:
C
C       * RES      : result of the integration of F over the whole 
C                       interval [1,N_END]
C
C
C  References: M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
C                "Functions", 9th Dover printing, pp.886-887, Dover 
C
C              P. A. Almeida Magalhaes Jr and C. Almeida Magalhaes,
C                J. Math. Stat. 6, 193-204 (2010)
C
C  This version: closed Newton-Cotes formula limited to N_RULE = 6
C                no open Newton-Cotes formula included
C                Euler-MacLaurin formula limited to N_RULE = 3
C                                                         
C   Author :  D. Sébilleau
C
C                                          Last modified : 29 Oct 2014
C     
C
      REAL*4 F(N_POINTS),F_1(N_POINTS),F_3(N_POINTS),F_5(N_POINTS)
      REAL*4 X(N_POINTS)
      REAL*4 F_INT1,F_INT2,F_INT3,F_INT4
      REAL*4 RES,RES0,RES1,C_H
      REAL*4 BN(0:6),H,H1
      REAL*4 CNC2(2),CNC3(3),CNC4(4),CNC5(5),CNC6(6)
      REAL*4 CN(6)
      REAL*4 P,A(10)
C
      CHARACTER*3 METH
C
C  Bernouilli numbers
C
      DATA BN /1.,-0.5,0.166667,0.,0.033333,0.,0.023809/
C
C  Closed formula Newton-Cotes coefficients CNCn for n-point formula
C
      DATA CNC2  /1.0,1.0/
      DATA CNC3  /1.0,4.0,1.0/
      DATA CNC4  /1.0,3.0,3.0,1.0/
      DATA CNC5  /7.0,32.0,12.0,32.0,7.0/
      DATA CNC6  /19.0,75.0,50.0,50.0,75.0,19.0/
C
      DATA CN    /0.0,2.0,6.0,8.0,90.0,288.0/
C
C  Checking for consitency of input data
C
      IF(N_BEG.LT.1) THEN
        WRITE(6,10)
        STOP
      ENDIF
      IF(N_END.GT.N_POINTS) THEN
        WRITE(6,20)
        STOP
      ENDIF
      IF(METH.EQ.'NCQ') THEN
        IF((N_RULE.LT.2).OR.(N_RULE.GT.6)) THEN
          WRITE(6,30)
          STOP
        ENDIF
      ELSEIF(METH.EQ.'EMS') THEN
        IF((N_RULE.LT.2).OR.(N_RULE.GT.4)) THEN
          WRITE(6,40)
          STOP
        ENDIF
      ENDIF        
C
      H=X(2)-X(1)
C
      I_FLAG=N_RULE-1
C
C
C  Computation of Int_{1}^{X} F(X) dX for X in [N_BEG,N_END]
C
C
C  The number of points used for each 
C    formula is N_RULE. (N_END-N_BEG-1) must
C    must be divisible by I_FLAG in 
C    order to fully apply the formula.
C    So, the formula is applied in 
C    the interval [N_BEG,N_END-N_REM],
C    where N_REM is the remainder of 
C    the division of (N_END-N_BEG-1) by I_FLAG,
C    and for the remaining interval,
C    an interpolation is used to 
C    obtain exactly I_FLAG+1 points
C    (F_INT1,F_INT2,F_INT3,F_INT4).
C    We note N_END-N_REM-1 = N_FIN.
C
      IF(METH.EQ.'NCQ') THEN
C
        N_REM=MOD(N_END-N_BEG,I_FLAG)
        N_FIN=N_END-N_REM-1
        C_H=FLOAT(I_FLAG)/CN(N_RULE)
        RES0=0.0
C
        IF(I_FLAG.EQ.1) THEN
C
C.............  2-point formula  ........
C 
          DO J=N_BEG,N_FIN,I_FLAG
            RES0=RES0+CNC2(1)*F(J)+CNC2(2)*F(J+1)
          ENDDO
          RES=RES0*H*C_H
C
        ELSEIF(I_FLAG.EQ.2) THEN
C
C.............  3-point formula  ........
C 
          IF(N_FIN.GT.N_BEG) THEN
            DO J=N_BEG,N_FIN,I_FLAG
              RES0=RES0+CNC3(1)*F(J)+CNC3(2)*F(J+1)+CNC3(3)*F(J+2)
            ENDDO
          ENDIF
          RES0=RES0*H*C_H
C
          IF(N_REM.EQ.0) THEN
            RES=RES0
          ELSEIF(N_REM.EQ.1) THEN
C
C               Lagrange 3-point interpolation for step H/2 point
C                  (or Lagrange 2-point when not possible)
C
            P=1.0/2.0
            IF(N_END.GT.2) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT1=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT1=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            RES1=CNC3(1)*F(N_END-1)+CNC3(2)*F_INT1+CNC3(3)*F(N_END)
            H1=H/2.0
            RES=RES0+RES1*H1*C_H
          ENDIF
C
        ELSEIF(I_FLAG.EQ.3) THEN
C
C.............  4-point formula  ........
C 
          IF(N_FIN.GT.N_BEG) THEN
            DO J=N_BEG,N_FIN,I_FLAG
              RES0=RES0+CNC4(1)*F(J)+CNC4(2)*F(J+1)+CNC4(3)*F(J+2)+
     1                  CNC4(4)*F(J+3)
            ENDDO
          ENDIF
          RES0=RES0*H*C_H
C
          IF(N_REM.EQ.0) THEN
            RES=RES0
          ELSEIF(N_REM.EQ.1) THEN
C
C               Lagrange 3-point interpolation for step H/3 points
C                  (or Lagrange 2-point when not possible)
C
            P=1.0/3.0
            IF(N_END.GT.2) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT1=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT1=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            P=2.0/3.0
            IF(N_END.GT.2) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT2=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT2=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            RES1=CNC4(1)*F(N_END-1)+CNC4(2)*F_INT1+
     1           CNC4(3)*F_INT2+CNC4(4)*F(N_END)
            H1=H/3.0
            RES=RES0+RES1*H1*C_H
C
          ELSEIF(N_REM.EQ.2) THEN
C
C               Lagrange 3-point interpolation for step 2H/3 points
C                  (or Lagrange 2-point when not possible)
C               (F(N_END-1) is not used for the calculation of integral)
C
            P=2.0/3.0
            IF(N_END.GT.3) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT1=A(1)*F(N_END-3)+A(2)*F(N_END-2)+A(3)*F(N_END-1)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT1=A(1)*F(N_END-2)+A(2)*F(N_END-1)
            ENDIF
C
            P=1.0/3.0
            IF(N_END.GT.3) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT2=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT2=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            RES1=CNC4(1)*F(N_END-2)+CNC4(2)*F_INT1+CNC4(3)*F_INT2+
     1           CNC4(4)*F(N_END)
            H1=2.0*H/3.0
            RES=RES0+RES1*H1*C_H
C
          ENDIF
C
        ELSEIF(I_FLAG.EQ.4) THEN
C
C.............  5-point formula  ........
C 
          IF(N_FIN.GT.N_BEG) THEN
            DO J=N_BEG,N_FIN,I_FLAG
              RES0=RES0+CNC5(1)*F(J)+CNC5(2)*F(J+1)+CNC5(3)*F(J+2)+
     1                  CNC5(4)*F(J+3)+CNC5(5)*F(J+4)
            ENDDO
          ENDIF
          RES0=RES0*H*C_H
C
          IF(N_REM.EQ.0) THEN
            RES=RES0
          ELSEIF(N_REM.EQ.1) THEN
C
C               Lagrange 3-point interpolation for step H/4 points
C                  (or Lagrange 2-point when not possible)
C
            P=1.0/4.0
            IF(N_END.GT.2) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT1=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT1=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            P=1.0/2.0
            IF(N_END.GT.2) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT2=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT2=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            P=3.0/4.0
            IF(N_END.GT.2) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT3=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT3=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            RES1=CNC5(1)*F(N_END-1)+CNC5(2)*F_INT1+CNC5(3)*F_INT2+
     1           CNC5(4)*F_INT3+CNC5(5)*F(N_END)
            H1=H/4.0
            RES=RES0+RES1*H1*C_H
C
          ELSEIF(N_REM.EQ.2) THEN
C
C               Lagrange 3 point interpolation for step 2H/4 points
C                  (or Lagrange 2-point when not possible)
C
            P=1.0/2.0
            IF(N_END.GT.3) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT1=A(1)*F(N_END-3)+A(2)*F(N_END-2)+A(3)*F(N_END-1)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT1=A(1)*F(N_END-2)+A(2)*F(N_END-1)
            ENDIF
C
            P=1.0/2.0
            IF(N_END.GT.3) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT3=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT3=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            RES1=CNC5(1)*F(N_END-2)+CNC5(2)*F_INT1+
     1           CNC5(3)*F(N_END-1)+CNC5(4)*F_INT3+CNC5(5)*F(N_END)
            H1=H/2.0
            RES=RES0+RES1*H1*C_H
C
          ELSEIF(N_REM.EQ.3) THEN
C
C               Lagrange 3 point interpolation for step 3H/4 points
C
            P=3.0/4.0
            CALL LAGR_INTERP(3,P,A)
            F_INT1=A(1)*F(N_END-4)+A(2)*F(N_END-3)+A(3)*F(N_END-2)
C
            P=1.0/2.0
            CALL LAGR_INTERP(3,P,A)
            F_INT2=A(1)*F(N_END-3)+A(2)*F(N_END-2)+A(3)*F(N_END-1)
C
            P=1.0/4.0
            CALL LAGR_INTERP(3,P,A)
            F_INT3=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
C
            RES1=CNC5(1)*F(N_END-3)+CNC5(2)*F_INT1+
     1           CNC5(3)*F_INT2+CNC5(4)*F_INT3+CNC5(5)*F(N_END)
            H1=3.0*H/4.0
            RES=RES0+RES1*H1*C_H
C
          ENDIF
C
        ELSEIF(I_FLAG.EQ.5) THEN
C
C.............  6-point formula  ........
C 
          IF(N_FIN.GT.N_BEG) THEN
            DO J=N_BEG,N_FIN,I_FLAG
              RES0=RES0+CNC6(1)*F(J)+CNC6(2)*F(J+1)+CNC6(3)*F(J+2)+
     1                  CNC6(4)*F(J+3)+CNC6(5)*F(J+4)+CNC6(6)*F(J+5)
            ENDDO
          ENDIF
          RES0=RES0*H*C_H
C
          IF(N_REM.EQ.0) THEN
            RES=RES0
          ELSEIF(N_REM.EQ.1) THEN
C
C               Lagrange 3-point interpolation for step H/5 points
C                  (or Lagrange 2-point when not possible)
C
            P=1.0/5.0
            IF(N_END.GT.2) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT1=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT1=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            P=2.0/5.0
            IF(N_END.GT.2) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT2=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT2=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            P=3.0/5.0
            IF(N_END.GT.2) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT3=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT3=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            P=4.0/5.0
            IF(N_END.GT.2) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT4=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT4=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            RES1=CNC6(1)*F(N_END-1)+CNC6(2)*F_INT1+CNC6(3)*F_INT2+
     1           CNC6(4)*F_INT3+CNC6(5)*F_INT4+CNC6(6)*F(N_END)
            H1=H/5.0
            RES=RES0+RES1*H1*C_H
C
          ELSEIF(N_REM.EQ.2) THEN
C
C               Lagrange 3 point interpolation for step 2H/5 points
C                  (or Lagrange 2-point when not possible)
C
            P=2.0/5.0
            IF(N_END.GT.3) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT1=A(1)*F(N_END-3)+A(2)*F(N_END-2)+A(3)*F(N_END-1)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT1=A(1)*F(N_END-2)+A(2)*F(N_END-1)
            ENDIF
C
            P=4.0/5.0
            IF(N_END.GT.3) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT2=A(1)*F(N_END-3)+A(2)*F(N_END-2)+A(3)*F(N_END-1)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT2=A(1)*F(N_END-2)+A(2)*F(N_END-1)
            ENDIF
C
            P=1.0/5.0
            IF(N_END.GT.3) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT3=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT3=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            P=3.0/5.0
            IF(N_END.GT.3) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT4=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT4=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            RES1=CNC6(1)*F(N_END-2)+CNC6(2)*F_INT1+CNC6(3)*F_INT2+
     1           CNC6(4)*F_INT3+CNC6(5)*F_INT4+CNC6(6)*F(N_END)
            H1=2.0*H/5.0
            RES=RES0+RES1*H1*C_H
C
          ELSEIF(N_REM.EQ.3) THEN
C
C               Lagrange 3 point interpolation for step 3H/5 points
C                  (or Lagrange 2-point when not possible)
C
            P=3.0/5.0
            IF(N_END.GT.4) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT1=A(1)*F(N_END-4)+A(2)*F(N_END-3)+A(3)*F(N_END-2)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT1=A(1)*F(N_END-3)+A(2)*F(N_END-2)
            ENDIF
C
            P=1.0/5.0
            IF(N_END.GT.4) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT2=A(1)*F(N_END-3)+A(2)*F(N_END-2)+A(3)*F(N_END-1)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT2=A(1)*F(N_END-2)+A(2)*F(N_END-1)
            ENDIF
C
            P=4.0/5.0
            IF(N_END.GT.4) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT3=A(1)*F(N_END-3)+A(2)*F(N_END-2)+A(3)*F(N_END-1)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT3=A(1)*F(N_END-2)+A(2)*F(N_END-1)
            ENDIF
C
            P=2.0/5.0
            IF(N_END.GT.4) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT4=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT4=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            RES1=CNC6(1)*F(N_END-3)+CNC6(2)*F_INT1+CNC6(3)*F_INT2+
     1           CNC6(4)*F_INT3+CNC6(5)*F_INT4+CNC6(6)*F(N_END)
            H1=3.0*H/5.0
            RES=RES0+RES1*H1*C_H
C
          ELSEIF(N_REM.EQ.4) THEN
C
C               Lagrange 3 point interpolation for step 4H/5 points
C                  (or Lagrange 2-point when not possible)
C
            P=4.0/5.0
            IF(N_END.GT.5) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT1=A(1)*F(N_END-5)+A(2)*F(N_END-4)+A(3)*F(N_END-3)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT1=A(1)*F(N_END-4)+A(2)*F(N_END-3)
            ENDIF
C
            P=3.0/5.0
            IF(N_END.GT.5) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT2=A(1)*F(N_END-4)+A(2)*F(N_END-3)+A(3)*F(N_END-2)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT2=A(1)*F(N_END-3)+A(2)*F(N_END-2)
            ENDIF
C
            P=2.0/5.0
            IF(N_END.GT.5) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT3=A(1)*F(N_END-3)+A(2)*F(N_END-2)+A(3)*F(N_END-1)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT3=A(1)*F(N_END-2)+A(2)*F(N_END-1)
            ENDIF
C
            P=1.0/5.0
            IF(N_END.GT.5) THEN
              CALL LAGR_INTERP(3,P,A)
              F_INT4=A(1)*F(N_END-2)+A(2)*F(N_END-1)+A(3)*F(N_END)
            ELSE
              CALL LAGR_INTERP(2,P,A)
              F_INT4=A(1)*F(N_END-1)+A(2)*F(N_END)
            ENDIF
C
            RES1=CNC6(1)*F(N_END-4)+CNC6(2)*F_INT1+CNC6(3)*F_INT2+
     1           CNC6(4)*F_INT3+CNC6(5)*F_INT4+CNC6(6)*F(N_END)
            H1=4.0*H/5.0
            RES=RES0+RES1*H1*C_H
C
            ENDIF
C
        ENDIF
C
      ELSEIF(METH.EQ.'EMS') THEN
C
        IF(N_RULE.GE.1) THEN
          RES1=(F(N_BEG)+F(N_END))*0.5
          DO J=N_BEG+1,N_END-1
            RES1=RES1+F(J)
          ENDDO
          RES1=RES1*H
        ENDIF
        IF(N_RULE.GE.2) THEN
          RES1=RES1-BN(2)*H*H*(F_1(N_END)-F_1(N_BEG))/2.0
        ENDIF
        IF(N_RULE.GE.3) THEN
          RES1=RES1-BN(4)*H*H*H*H*(F_3(N_END)-F_3(N_BEG))/24.0
        ENDIF
        IF(N_RULE.GE.4) THEN
          RES1=RES1-BN(6)*H*H*H*H*H*H*(F_5(N_END)-F_5(N_BEG))/720.0
        ENDIF
        RES=RES1
C
      ENDIF
C
C  Formats
C
  10  FORMAT(//,10X,'<<<<<  WRONG VALUE OF N_BEG:   >>>>>',/,
     1          10X,'<<<<<  CANNOT BE LOWER THAN 1  >>>>>',//)
  20  FORMAT(//,10X,'<<<<<  WRONG VALUE OF N_END:   >>>>>',/,
     1          10X,'<<<<<  CANNOT EXCEED N_POINTS  >>>>>',//)
  30  FORMAT(//,10X,'<<<<<  WRONG VALUE OF N_RULE:  >>>>>',/,
     1          10X,'<<<<<    SHOULD BE IN [2,6]    >>>>>',//)
  40  FORMAT(//,10X,'<<<<<  WRONG VALUE OF N_RULE:  >>>>>',/,
     1          10X,'<<<<<    SHOULD BE IN [2,4]    >>>>>',//)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE LAGR_INTERP(N,P,A)
C
C      This subroutine computes the coefficients for the Lagrange 
C        n-point interpolation, 1 < n < 7
C
C  Input parameters:
C
C       * N        : number of points of the interpolation
C       * P        : value of the step fraction
C
C  Output parameters:
C
C       * A(N)     : coefficients of the interpolation
C
C
C  References: M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
C                "Functions", 9th Dover printing, pp.878-879, Dover 
C
C                                                         
C   Author :  D. Sébilleau
C
C                                          Last modified : 29 Oct 2014
C   
      INTEGER N
C
      REAL*4 P,A(10)
C
C  Initialization
C
      DO J=1,10
         A(J)=0.
      ENDDO
C
      IF(N.EQ.2) THEN
C
C..........  2-point Lagrange interpolation  ............
C
        A(1)=1.0-P
        A(2)=P
C
      ELSEIF(N.EQ.3) THEN
C
C..........  3-point Lagrange interpolation  ............
C
        A(1)=0.5*P*(P-1.0)
        A(2)=1.0-P*P
        A(3)=0.5*P*(P+1.0)
C
      ELSEIF(N.EQ.4) THEN
C
C..........  4-point Lagrange interpolation  ............
C
        A(1)=-P*(P-1.0)*(P-2.0)/6.0
        A(2)=(P*P-1.0)*(P-2.0)/2.0
        A(3)=-P*(P+1.0)*(P-2.0)/2.0
        A(4)=P*(P*P-1.0)/6.0
C
      ELSEIF(N.EQ.5) THEN
C
C..........  5-point Lagrange interpolation  ............
C
        A(1)=(P*P-1.0)*P*(P-2.0)/24.0
        A(2)=-(P-1.0)*P*(P*P-4.0)/6.0
        A(3)=(P*P-1.0)*(P*P-4.0)/4.0
        A(4)=-(P+1.0)*P*(P*P-4.0)/6.0
        A(5)=(P*P-1.0)*P*(P+2.0)/24.0
C
      ELSEIF(N.EQ.6) THEN
C
C..........  6-point Lagrange interpolation  ............
C
        A(1)=-P*(P*P-1.0)*(P-2.0)*(P-3.0)/120.0
        A(2)=P*(P-1.0)*(P*P-4.0)*(P-3.0)/24.0
        A(3)=-(P*P-1.0)*(P*P-4.0)*(P-3.0)/12.0
        A(4)=P*(P+1.0)*(P*P-4.0)*(P-3.0)/12.0
        A(5)=-P*(P*P-1.0)*(P+2.0)*(P-3.0)/24.0
        A(6)=P*(P*P-1.0)*(P*P-4.0)/120.0
C
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE CHORDS(X,Y,N_POINTS,I_CHORD,METHOD,VALUE,N_BIN,IUO1,CH)
C
C  This subroutine computes the chords of a curve given 
C    by the arrays (X(I),Y(I)) with I = 1, N_POINTS
C
C
C  Input parameters:
C
C       X        :  x coordinates of the input file
C       Y        :  y coordinates of the input file
C       N_POINTS :  number of points in the file
C       I_CHORD  :  flag to select the type of chord calculation
C
C              = 1     chord_length from point I
C              = 2     distance from point I to chord (I-K)-(I+K)
C              = 3     chord length along direction THETA
C       METHOD   :  method used for the chords
C
C              = SIN   single value (given by VALUE)
C              = HIS   histogram
C              = SUM   sum/number of values
C       VALUE    :  value where the chord is to be computed
C
C             = n :
C                I_CHORD = 1 : all chords computed from point I = n
C                I_CHORD = 2 : for each point I, distance computed 
C                                from point I to chord Y(I-K)-Y(I+K)
C                I_CHORD = 3 : theta angle equal to n PI/32 
C                                with 0 < n < 16 
C       N_BIN    :  number of bins for the histogram (METHOD=HIS)
C
C
C  Output parameters:
C
C       CH       :  chord result
C
C
C   Author :  D. Sébilleau
C
C                                          Last modified : 11 Aug 2014
C
      PARAMETER (N_SIZE=1000)
C
      INTEGER VALUE
C
      REAL*4 X(N_SIZE),Y(N_SIZE),CH(N_SIZE),CHORD(N_SIZE,N_SIZE)
      REAL*4 X_I,Y_I,DIST,SUMJ,CHORD_MIN,CHORD_MAX,SIZE_BIN
      REAL*4 D_PLUS,A,SUMK
C
      CHARACTER*3 METHOD
C
      DATA PI /3.141593/
C
      X_POINTS=FLOAT(N_POINTS)
      N_MAX=MAX(N_POINTS,N_BIN)
C
C  Dimensionality checks
C
      IF(N_POINTS.GE.N_SIZE) THEN
        WRITE(IUO1,10) N_POINTS+1
        STOP
      ENDIF
C
      IF(METHOD.EQ.'HIS') THEN
        IF(N_BIN.GT.N_SIZE) THEN
          WRITE(IUO1,10) N_BIN
          STOP
        ENDIF
      ENDIF
C
      IF(I_CHORD.EQ.3) THEN
        IF((VALUE.LE.0).OR.(VALUE.GE.16)) THEN
          WRITE(IUO1,20)
        ENDIF
      ENDIF
C
      IF(METHOD.EQ.'SIN') THEN
        IF(VALUE.EQ.0) THEN
            WRITE(IUO1,30)
        ENDIF
      ENDIF
C
C  Initializations
C
      DO J=1,N_MAX
        CH(J)=0.0
      ENDDO
C
      DO I=1,N_POINTS
        DO J=1,N_POINTS
          CHORD(I,J)=0.0
        ENDDO
      ENDDO
C
C  Chord length from point I = VALUE
C
      IF(I_CHORD.EQ.1) THEN
C
        IF(METHOD.EQ.'SIN') THEN
C
          X_I=X(VALUE)
          Y_I=Y(VALUE)
          CH(VALUE)=0.0
          DO J=1,N_POINTS
            IF(J.EQ.VALUE) GOTO 5
            DIST=SQRT((X(J)-X_I)*(X(J)-X_I)+(Y(J)-Y_I)*(Y(J)-Y_I))
            CH(J)=DIST
    5       CONTINUE
          ENDDO
C
        ELSEIF(METHOD.EQ.'HIS')  THEN
C
C......... Computing chord lengths and the min/max
C
          CHORD_MIN=1.0E+30
          CHORD_MAX=0.00
          DO I=1,N_POINTS
            X_I=X(I)
            Y_I=Y(I)
            CHORD(I,I)=0.0
            SUMJ=0.0
            DO J=1,N_POINTS
              IF(J.EQ.I) GOTO 15
              DIST=SQRT((X(J)-X_I)*(X(J)-X_I)+(Y(J)-Y_I)*(Y(J)-Y_I))
              CHORD_MIN=MIN(DIST,CHORD_MIN)
              CHORD_MAX=MAX(DIST,CHORD_MAX)
              CHORD(I,J)=DIST
   15         CONTINUE
            ENDDO
          ENDDO
C
C......... Putting chords into bins
C
          SIZE_BIN=(CHORD_MAX-CHORD_MIN)/FLOAT(N_BIN)
C
          DO I=1,N_POINTS
            DO J=1,N_POINTS
              IF(J.EQ.I) GOTO 25
              DO L=1,N_BIN
                DIS_BIN_LO=CHORD_MIN+FLOAT(L-1)*SIZE_BIN
                DIS_BIN_UP=DIS_BIN_LO+SIZE_BIN
                JBIN=0
                IF(CHORD(I,J).GE.DIS_BIN_LO) JBIN=JBIN+1
                IF(CHORD(I,J).LE.DIS_BIN_UP) JBIN=JBIN+1
                IF(JBIN.EQ.2) THEN
                  CH(L)=CH(L)+1.0
                  GOTO 25
                ENDIF
              ENDDO
  25          CONTINUE
            ENDDO
          ENDDO
C
        ELSEIF(METHOD.EQ.'SUM') THEN
C
          DO I=1,N_POINTS
            X_I=X(I)
            Y_I=Y(I)
            SUMJ=0.0
            DO J=1,N_POINTS
              IF(J.EQ.I) GOTO 35
              DIST=SQRT((X(J)-X_I)*(X(J)-X_I)+(Y(J)-Y_I)*(Y(J)-Y_I))
              SUMJ=SUMJ+DIST
  35          CONTINUE
            ENDDO
            CH(I)=SUMJ/X_POINTS
          ENDDO
C  
        ENDIF
C
C  Distance from point I to chord (I-VALUE)-(I+VALUE)
C
      ELSEIF(I_CHORD.EQ.2) THEN
C
        IF(METHOD.EQ.'SIN') THEN
C
          K=VALUE
          DO I=1,N_POINTS
            IF((I+K).LE.N_POINTS) THEN
              M=I+K
            ELSE
              M=I+K-N_POINTS
            ENDIF
            IF((I-K).GE.1) THEN
              N=I-K
            ELSE
              N=N_POINTS+I-K
            ENDIF
            A=(Y(M)-Y(N))/(X(M)-X(N))
            D_PLUS=ABS(A*(X(M)-X(I))+Y(I)-Y(M))/SQRT(A*A+1.0)
            CH(I)=D_PLUS
          ENDDO
C
        ELSEIF(METHOD.EQ.'HIS')  THEN
C
C......... Computing chord lengths and the min/max
C
          CHORD_MIN=1.0E+30
          CHORD_MAX=0.00
          DO I=1,N_POINTS
            DO K=1,N_POINTS-1
              IF((I+K).LE.N_POINTS) THEN
                M=I+K
              ELSE
                M=I+K-N_POINTS
              ENDIF
              IF((I-K).GE.1) THEN
                N=I-K
              ELSE
                N=N_POINTS+I-K
              ENDIF
              A=(Y(M)-Y(N))/(X(M)-X(N))
              D_PLUS=ABS(A*(X(M)-X(I))+Y(I)-Y(M))/SQRT(A*A+1.0)
              CHORD(I,K)=D_PLUS
              CHORD_MIN=MIN(D_PLUS,CHORD_MIN)
              CHORD_MAX=MAX(D_PLUS,CHORD_MAX)
            ENDDO
          ENDDO
C
C......... Putting chords into bins
C
          SIZE_BIN=(CHORD_MAX-CHORD_MIN)/FLOAT(N_BIN)
C
          DO I=1,N_POINTS
            DO J=1,N_POINTS
              IF(J.EQ.I) GOTO 55
              DO L=1,N_BIN
                DIS_BIN_LO=CHORD_MIN+FLOAT(L-1)*SIZE_BIN
                DIS_BIN_UP=DIS_BIN_LO+SIZE_BIN
                JBIN=0
                IF(CHORD(I,J).GE.DIS_BIN_LO) JBIN=JBIN+1
                IF(CHORD(I,J).LE.DIS_BIN_UP) JBIN=JBIN+1
                IF(JBIN.EQ.2) THEN
                  CH(L)=CH(L)+1.0
                  GOTO 55
                ENDIF
              ENDDO
  55          CONTINUE
            ENDDO
          ENDDO
C
        ELSEIF(METHOD.EQ.'SUM') THEN
C
          NK=(N_POINTS-1)/2
          DO I=1,N_POINTS
            SUMK=0.0
            DO K=1,NK
              IF((I+K).LE.N_POINTS) THEN
                M=I+K
              ELSE
                M=I+K-N_POINTS
              ENDIF
              IF((I-K).GE.1) THEN
                N=I-K
              ELSE
                N=N_POINTS+I-K
              ENDIF
              A=(Y(M)-Y(N))/(X(M)-X(N))
              D_PLUS=ABS(A*(X(M)-X(I))+Y(I)-Y(M))/SQRT(A*A+1.0)
              SUMK=SUMK+D_PLUS
            ENDDO
            CH(I)=SUMK/FLOAT(NK)
          ENDDO
C
        ENDIF
C
C  Distance from point I to x axis along direction THETA = VALUE*PI/32
C
      ELSEIF(I_CHORD.EQ.3) THEN
C
        IF(METHOD.EQ.'SIN') THEN
C
          THETA=FLOAT(VALUE)*PI/32.
          SIN_TH=SIN(THETA)
          DO I=1,N_POINTS
            CH(I)=Y(I)*SIN_TH
          ENDDO
C
        ELSEIF(METHOD.EQ.'HIS') THEN
C
          CHORD_MIN=1.0E+30
          CHORD_MAX=0.00
          DO I=1,N_POINTS
            DO K=1,15
              THETA=FLOAT(K)*PI/32.
              SIN_TH=SIN(THETA)
              DIST=Y(I)*SIN_TH
              CHORD(I,K)=DIST
              CHORD_MIN=MIN(DIST,CHORD_MIN)
              CHORD_MAX=MAX(DIST,CHORD_MAX)
            ENDDO
          ENDDO
C
C......... Putting chords into bins
C
          SIZE_BIN=(CHORD_MAX-CHORD_MIN)/FLOAT(N_BIN)
C
          DO I=1,N_POINTS
            DO J=1,15
              IF(J.EQ.I) GOTO 75
              DO L=1,N_BIN
                DIS_BIN_LO=CHORD_MIN+FLOAT(L-1)*SIZE_BIN
                DIS_BIN_UP=DIS_BIN_LO+SIZE_BIN
                JBIN=0
                IF(CHORD(I,J).GE.DIS_BIN_LO) JBIN=JBIN+1
                IF(CHORD(I,J).LE.DIS_BIN_UP) JBIN=JBIN+1
                IF(JBIN.EQ.2) THEN
                  CH(L)=CH(L)+1.0
                  GOTO 75
                ENDIF
              ENDDO
  75          CONTINUE
            ENDDO
          ENDDO
C
        ELSEIF(METHOD.EQ.'SUM') THEN
C
          DO I=1,N_POINTS
            SUMK=0.0
            DO K=1,15
              THETA=FLOAT(K)*PI/32.
              SIN_TH=SIN(THETA)
              SUMK=SUMK+Y(I)*SIN_TH
            ENDDO
            CH(I)=SUMK/15.0
          ENDDO
C
        ENDIF
C
      ENDIF
C
C  Writing the transformed coordinates
C
      IF(METHOD.EQ.'HIS')  THEN
        DO J=1,N_BIN
          WRITE(98,*) J,CH(J)
        ENDDO
      ELSE
        DO J=1,N_POINTS
          WRITE(98,*) J,CH(J)
        ENDDO
      ENDIF
C
C  Formats
C
  10  FORMAT(//,10X,'<<<<<  N_SIZE SHOULD BE LARGER THAN ',I5,
     1      '  >>>>>',//)
  20  FORMAT(//,10X,'<<<<<  ERROR IN THE INPUT DATA FILE:  >>>>>',/
     1              '<<<<<    VALUE SHOULD BE IN ]0,16[    >>>>>',//)
  30  FORMAT(//,10X,'<<<<<  ERROR IN THE INPUT DATA FILE:  >>>>>',/
     1              '<<<<<      VALUE CANNOT BE ZERO       >>>>>',//)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE CHAIN_CODE(X,Y,N_CMP,N_CONNECT,SCALEC,IUO1,CC)
C
C  This subroutine computes the chain code of a curve given 
C    by the arrays (X(I),Y(I))
C
C
C  Input parameters:
C
C       X           :  x coordinates of the input file
C       Y           :  y coordinates of the input file
C       N_CMP       :  number of points in the file
C       N_CONNECT   :  connectivity of the chain code (can be 3, 5 or 9)
C       SCALEC      :  scaling factor to compute the tangent angle
C       IUO1        :  checkfile index for printing
C
C
C  Output parameters:
C
C       CC          :  chain code
C
C
C   Author :  D. Sébilleau
C
C                                          Last modified : 10 Sep 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 Y(N_SIZE),X(N_SIZE),TG
      REAL*4 X_MIN,X_MAX,Y_MIN,Y_MAX,SCALEC,SCALE1
C
      INTEGER CC(N_SIZE)
C
      DATA TAN_PIO4,TAN_PIO8,TAN_3PIO8       / 1.0,0.414214,2.414214/
      DATA TAN_MPIO4,TAN_MPIO8,TAN_M3PIO8    /-1.0,-0.414214,-2.414214/
      DATA TAN_PIO16,TAN_3PIO16              /0.198912,0.668179/
      DATA TAN_5PIO16,TAN_7PIO16             /1.496606,5.027339/
      DATA TAN_MPIO16,TAN_M3PIO16            /-0.198912,-0.668179/
      DATA TAN_M5PIO16,TAN_M7PIO16           /-1.496606,-5.027339/
C
C  Dimensionality check
C
      IF((N_CMP+1).GE.N_SIZE) THEN
        WRITE(IUO1,10) N_CMP+2
        STOP
      ENDIF
C
C  Connectivity check
C
      I_CHK=0
      IF(N_CONNECT.EQ.3) THEN
        I_CHK=I_CHK+1
      ELSEIF(N_CONNECT.EQ.5) THEN
        I_CHK=I_CHK+1
      ELSEIF(N_CONNECT.EQ.9) THEN
        I_CHK=I_CHK+1
      ENDIF
C
      IF(I_CHK.EQ.0) THEN
        WRITE(IUO1,20) 
        STOP
      ENDIF
C
C  Setting point (N_CMP+1) = 1 to have 
C      a N_CMP long chain
C
      X(N_CMP+1)=X(1)
      Y(N_CMP+1)=Y(1)  
C
C  Computing the automatic scaling factor used to have 
C             meaningful tangent angles 
C
      X_MIN=1.E+30
      Y_MIN=1.E+30
      X_MAX=-1.E+30
      Y_MAX=-1.E+30
C
      DO I=1,N_CMP
C
        X_MIN=MIN(X_MIN,X(I))
        X_MAX=MAX(X_MAX,X(I))
        Y_MIN=MIN(Y_MIN,Y(I))
        Y_MAX=MAX(Y_MAX,Y(I))
C
      ENDDO
C
      SCALE1=(X_MAX-X_MIN)/(Y_MAX-Y_MIN)
      SCALE1=SCALE1*SCALEC
      WRITE(IUO1,30) SCALE1
C
C  Computation of the chain code
C
      DO I=1,N_CMP
C
        TG=SCALE1*(Y(I+1)-Y(I))/(X(I+1)-X(I))
C
        IF(N_CONNECT.EQ.3) THEN
C
C  3-connectivity clockwise chain code
C
          IF(TG.GE.TAN_PIO4) THEN
            CC(I)=0
          ELSEIF(TG.LE.TAN_MPIO4) THEN
            CC(I)=2
          ELSE
            CC(I)=1
          ENDIF
C
        ELSEIF(N_CONNECT.EQ.5) THEN
C
C  5-connectivity clockwise chain code
C    (dividing into 2 quadrants)
C
          IF(TG.GE.0.0) THEN
C
C.............. Top quadrant [0,PI/2[
C
            IF(TG.GE.TAN_3PIO8) THEN
              CC(I)=0
            ELSEIF(TG.LE.TAN_PIO8) THEN
              CC(I)=2
            ELSE
              CC(I)=1
            ENDIF
          ELSE
C
C.............. Bottom quadrant [0,-PI/2[
C
            IF(TG.GE.TAN_MPIO8) THEN
              CC(I)=2
            ELSEIF(TG.LE.TAN_M3PIO8) THEN
              CC(I)=4
            ELSE
              CC(I)=3
            ENDIF
C
          ENDIF         
C
        ELSEIF(N_CONNECT.EQ.9) THEN
C
C  9-connectivity clockwise chain code
C
          IF(TG.GE.TAN_7PIO16) THEN
            CC(I)=0
          ELSEIF((TG.LT.TAN_7PIO16).AND.(TG.GE.TAN_5PIO16)) THEN
            CC(I)=1
          ELSEIF((TG.LT.TAN_5PIO16).AND.(TG.GE.TAN_3PIO16)) THEN
            CC(I)=2
          ELSEIF((TG.LT.TAN_3PIO16).AND.(TG.GE.TAN_PIO16)) THEN
            CC(I)=3
          ELSEIF((TG.LT.TAN_PIO16).AND.(TG.GE.TAN_MPIO16)) THEN
            CC(I)=4
          ELSEIF((TG.LT.TAN_MPIO16).AND.(TG.GE.TAN_M3PIO16)) THEN
            CC(I)=5
          ELSEIF((TG.LT.TAN_M3PIO16).AND.(TG.GE.TAN_M5PIO16)) THEN
            CC(I)=6
          ELSEIF((TG.LT.TAN_M5PIO16).AND.(TG.GE.TAN_M7PIO16)) THEN
            CC(I)=7
          ELSE
            CC(I)=8
          ENDIF 
C        
        ENDIF
C
      ENDDO
C
C  Writing the chain code into a file
C
      DO I=1,N_CMP
        WRITE(98,*) I,CC(I)
      ENDDO
C
C  Formats
C
  10  FORMAT(//,10X,'<<<<<  N_SIZE SHOULD BE LARGER THAN ',I5,
     1      '  >>>>>',//)
  20  FORMAT(//,10X,'<<<<<  N_CONNECT SHOULD BE 3, 5 or 9  >>>>>',//)
  30  FORMAT(6X,'--->  A SCALING FACTOR OF ',E12.6,' IS APPLIED ',
     1           'TO HAVE A MEANINGFUL TANGENT ANGLE')
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE CONTOUR(I,X,N_POINTS,IUO1,N_BIN,N_LEN,SH_AN,I_FOU,
     1                   I_METHOD,FILE,SA)
C
C  This subroutine transforms the points of a curve (I,X) into 
C    the shape descriptors SA described by SH_AN of a closed contour 
C    for further use
C
C
C  Input parameters:
C
C       I           :  y coordinates of the input file
C       X           :  x coordinates of the input file
C       N_POINTS    :  number of points in the file
C       IUO1        :  checkfile index for printing
C       N_BIN       :  number of bins (for histograms)
C       N_LEN       :  value at which the shape descriptor is computed
C                        when it is double-valued
C       SH_AN       :  type of shape descriptor used
C
C                    = CDIS  centroid distance
C                    = TANG  tangent angle
C                    = CURV  curvature function
C                    = TRAR  triangle area
C                    = BEAS  beam angle statistics
C                    = 8CCH  8-connectivity chain code
C                    = CLEN  chord length
C                    = CANG  chord angle
C                    = ACDI  arc chord distance
C                    = FOUR  Fourier descriptors
C       I_FOU       :  type of the real Fourier descriptor
C
C                    = 1 : modulus of the complex Fourier descriptor
C                    = 2 : argument of the complex Fourier descriptor 
C                    = 3 : real part of the complex Fourier descriptor
C                    = 4 : imaginary part of the complex Fourier descriptor
C       I_METHOD    :  normalization method for the contour
C
C                    = 1 : second order moments equal to unity 
C                    = 2 : affine orthogonalization method 
C                    = 3 : equal area method (area = 1) for contour
C                    = 4 : Gu-Kundu method
C       FILE        :  name of input file
C
C
C  Output parameters:
C
C       SA          :  shape descriptors
C
C
C   Author :  D. Sébilleau 
C
C                                          Last modified : 14 Aug 2014
C
      PARAMETER (N_SIZE=1000,NMAX=9999)
C
      REAL*4 I(N_SIZE),X(N_SIZE)
      REAL*4 XX(N_SIZE),YY(N_SIZE)
      REAL*4 SA(0:NMAX)
C
      CHARACTER*4 SH_AN
      CHARACTER*48 FILE
C
C  Transforming the curve into a closed contour
C
      CALL CURVE_TO_CONTOUR(I,X,N_POINTS,IUO1,I_METHOD,XX,YY)
C
C  Computing the shape descriptor SA
C
      CALL SHAPE_DESCRIPTORS(XX,YY,N_POINTS,N_BIN,N_LEN,IUO1,SH_AN,
     1                       FILE,I_FOU,SA)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE CURVE_TO_CONTOUR(I,X,N_POINTS,IUO1,I_METHOD,XX,YY)
C
C  This subroutine transforms a 2D curve I(N_POINTS),X(N_POINTS)
C    into a closed contour XX(N_POINTS),YY(N_POINTS)
C
C  The center of the contour is taken as (0,0) and the countour 
C    points ZZ(J) = XX(J) + i YY(J) are given by:
C
C          ZZ(J) = R(J) EXP(i TH(J))
C
C              with    R(J)  = Y(J)
C                      TH(J) = 2 PI * (X(J)-X(1))/(X(N_POINTS)-X(1) * 
C                                     N/(N+1)
C
C  Then, the contour is normalized to be scale-independent
C
C
C  Input parameters:
C
C       I           :  y coordinates of the input file
C       X           :  x coordinates of the input file
C       N_POINTS    :  number of points in the file
C       IUO1        :  checkfile index for printing
C       I_METHOD    :  normalization method for the contour
C
C                    = 1 : second order moments equal to unity 
C                    = 2 : affine orthogonalization method 
C                    = 3 : equal area method (area = 1) for contour
C                    = 4 : Gu-Kundu method
C
C
C  Output parameters:
C
C       XX          :  x coordinates of the contour
C       YY          :  y coordinates of the contour
C
C
C   Author :  D. Sébilleau 
C
C                                          Last modified : 14 Aug 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 I(N_SIZE),X(N_SIZE)
      REAL*4 XX(N_SIZE),YY(N_SIZE)
      REAL*4 ANGLE,RATIO
C
      DATA PI /3.141593/
C
      X_POINTS=FLOAT(N_POINTS)
C
C  Defining the contour
C
      DO J=1,N_POINTS
C
        RATIO=(X(J)-X(1))/(X(N_POINTS)-X(1))
        ANGLE=2.*PI*RATIO*X_POINTS/(X_POINTS+1.0)
        XX(J)=I(J)*COS(ANGLE)
        YY(J)=I(J)*SIN(ANGLE)
C
      ENDDO
C
C  Normalizing the contour whenever required
C
      IF(I_METHOD.EQ.0) THEN
        DO J=1,N_POINTS
          WRITE(99,*) XX(J),YY(J)
        ENDDO
        WRITE(99,*) XX(1),YY(1)
      ELSE
        CALL NORMALIZE_CONTOUR(XX,YY,N_POINTS,I_METHOD,IUO1)
        DO J=1,N_POINTS
          WRITE(99,*) XX(J),YY(J)
        ENDDO
        WRITE(99,*) XX(1),YY(1)
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE SHAPE_DESCRIPTORS(X,Y,N_POINTS,N_BIN,N_LEN,IUO1,
     1                             SH_AN,FILE,I_FOU,SA)
C
C  This routine computes various shape descriptors of a contour
C    containing N_POINTS (X,Y). By contour, we mean a closed curve
C    so that if the points vary from 1 to N_POINTS, we have 
C    X(N_POINTS+1) = X(1) and Y(N_POINTS+1) = Y(1). N_BIN is the 
C    number of bins for histograms (chord length and angle distributions)
C
C  References:  M. Yang, K. Kpalma and J. Ronsin, 
C               Pattern Recognition Techniques, Technology and Applications,
C               edited by Peng-Yen Yin, p. 626, I-Tech, Vienna (2008)
C
C
C  Input parameters:
C
C       I           :  y coordinates of the input file
C       X           :  x coordinates of the input file
C       N_POINTS    :  number of points in the file
C       N_BIN       :  number of bins (for histograms)
C       N_LEN       :  value at which the shape descriptor is computed
C                        when it is double-valued
C       IUO1        :  checkfile index for printing
C       SH_AN       :  type of shape descriptor used
C
C                    = CDIS  centroid distance
C                    = TANG  tangent angle
C                    = CURV  curvature function
C                    = TRAR  triangle area
C                    = BEAS  beam angle statistics
C                    = 8CCH  8-connectivity chain code
C                    = CLEN  chord length
C                    = CANG  chord angle
C                    = ACDI  arc chord distance
C                    = FOUR  Fourier descriptors
C       FILE        :  name of input file
C       I_FOU       :  type of the real Fourier descriptor
C
C                    = 1 : modulus of the complex Fourier descriptor
C                    = 2 : argument of the complex Fourier descriptor 
C                    = 3 : real part of the complex Fourier descriptor
C                    = 4 : imaginary part of the complex Fourier descriptor
C
C
C  Output parameters:
C
C       SA          :  shape descriptors
C
C
C   Author :  D. Sébilleau
C
C
C                                          Last modified : 14 Aug 2014
C
      PARAMETER (N_SIZE=1000,NMAX=9999)
C
      COMPLEX*8 Z(N_SIZE),AN(0:NMAX),ZEROC,IC,SUMC
C
      REAL*4 X(N_SIZE),Y(N_SIZE)
      REAL*4 D1X(N_SIZE),D2X(N_SIZE),D1Y(N_SIZE),D2Y(N_SIZE)
      REAL*4 F3(N_SIZE),F4(N_SIZE),F5(N_SIZE)
      REAL*4 NUM,DEN,TMP1,TMP2
      REAL*4 C_DIST(N_SIZE),T_ANGLE(N_SIZE),T2_ANGLE(N_SIZE)
      REAL*4 K(N_SIZE),K_NORM(N_SIZE)
      REAL*4 TAR(N_SIZE,N_SIZE)
      REAL*4 BAS(N_SIZE,N_SIZE)
      REAL*4 CHORD_DIS(N_SIZE,N_SIZE),CHORD_ANG(N_SIZE,N_SIZE),CD
      REAL*4 ACD(N_SIZE)
      REAL*4 SUM1,SUM2,SUM3,SUM4,A,B,C,ALPHA,D2IDA2
      REAL*4 PERIM,C_AREA,CENTROID(2)
      REAL*4 CXX,CXY,CYX,CYY,LAMBDA1,LAMBDA2,ECCEN,COVAR(2,2)
      REAL*4 CIRC_RATIO,DI,MUR,SIGMAR,COMPACT
      REAL*4 BE,TANGENT
      REAL*4 CHORD_DIS_MIN,CHORD_DIS_MAX,DIS_BIN_LO,DIS_BIN_UP
      REAL*4 SIZE_BIN,ANGL_BIN
      REAL*4 SA(0:NMAX)
      REAL*4 XI,XN,XAI,YAI,XX,YY
C
      INTEGER CHAIN_CODE(N_SIZE)
      INTEGER CHORD_DIS_HIS(N_SIZE),CHORD_ANG_HIS(N_SIZE)
C
      CHARACTER*4 SH_AN
      CHARACTER*48 FILE
      CHARACTER*50 FILE2
C
      DATA EPS                              /0.001/
      DATA TWOPI,PI,PIO2                    /6.283185,3.141593,1.570796/
      DATA TAN_PIO4,TAN_PIO8,TAN_3PIO8      / 1.0,0.414214,2.414214/
      DATA TAN_MPIO4,TAN_MPIO8,TAN_M3PIO8   /-1.0,-0.414214,-2.414214/
      DATA TAN_3PIO4,TAN_5PIO8,TAN_7PIO8    /-1.0,-2.414214,-0.414214/
      DATA TAN_M3PIO4,TAN_M5PIO8,TAN_M7PIO8 / 1.0,2.414214,0.414214/
C
      ZEROC=(0.,0.)
      IC=(0.0,1.0)
C
C  Real size of filename FILE
C
      N_CHAR=0
      DO J_CHAR=1,48
        IF(FILE(J_CHAR:J_CHAR).EQ.' ') GOTO 500
        N_CHAR=N_CHAR+1
      ENDDO
 500  CONTINUE
      FILE2=FILE(1:N_CHAR)//' :'
C
C  Dimensionality check
C
      IF(N_POINTS.GE.N_SIZE) THEN
        WRITE(IUO1,10) N_POINTS+1
        STOP
      ENDIF
      IF(SH_AN.EQ.'FOUR') THEN
        IF(N_BIN.GT.NMAX) THEN
          WRITE(IUO1,11) N_BIN+1
          STOP
        ENDIF
        IF(N_BIN.GT.N_POINTS) THEN
          N_BIN=N_POINTS
          WRITE(IUO1,13)
        ENDIF
      ENDIF
      IF(SH_AN.EQ.'ACDI') THEN
        IF(N_BIN.EQ.0) THEN
          WRITE(IUO1,12)
          STOP
        ENDIF
      ENDIF
C
C  Initializations
C
      DO J=1,N_BIN
        CHORD_DIS_HIS(J)=0
        CHORD_ANG_HIS(J)=0
      ENDDO
C
      DO J=0,NMAX
        SA(J)=0.
      ENDDO
C
      X_POINTS=FLOAT(N_POINTS)
C
C  Boundary condition for contour
C
      X(N_POINTS+1)=X(1)
      Y(N_POINTS+1)=Y(1)
C
C  Perimeter
C
      PERIM=0.0
      DO I=1,N_POINTS
         PERIM=PERIM+SQRT((X(I+1)-X(I))*(X(I+1)-X(I))+
     1                    (Y(I+1)-Y(I))*(Y(I+1)-Y(I)))
      ENDDO
C
C  Contour area C_AREA
C
      SUM1=0.0
      DO I=1,N_POINTS
         SUM1=SUM1+(X(I)*Y(I+1)-X(I+1)*Y(I))
      ENDDO
      C_AREA=0.5*ABS(SUM1)
C
C  Compactness
C
      COMPACT=2.*SQRT(C_AREA*PI)/PERIM
C
C  Centroid (CENTROID(1),CENTROID(2))
C
      SUM1=0.0
      SUM2=0.0
      DO I=1,N_POINTS
        SUM1=SUM1+(X(I)+X(I+1))*(X(I)*Y(I+1)-X(I+1)*Y(I))
        SUM2=SUM2+(Y(I)+Y(I+1))*(X(I)*Y(I+1)-X(I+1)*Y(I))
      ENDDO
      CENTROID(1)=SUM1/(6.*C_AREA)
      CENTROID(2)=SUM2/(6.*C_AREA)
C
C  Axis of least inertia angle slope angle THETA
C
      A=0.0
      B=0.0
      C=0.0
      DO I=1,N_POINTS
        A=A+X(I)*X(I)
        B=B+2.0*X(I)*Y(I)
        C=C+Y(I)*Y(I)
      ENDDO
      ALPHA=0.5*ATAN(B/(A-C))
      D2IDA2=2.0*(A-C)*COS(ALPHA+ALPHA)+2.0*B*SIN(ALPHA+ALPHA)
      IF(D2IDA2.LT.0.0) THEN
        THETA=(ALPHA+PIO2)*180.0/PI
      ELSE
        THETA=ALPHA*180.0/PI
      ENDIF
C
C  Covariance matrix COVAR and eccentricity ECCEN
C
      SUM1=0.0
      SUM2=0.0
      SUM3=0.0
      SUM4=0.0
      DO I=1,N_POINTS
        SUM1=SUM1+(X(I)-CENTROID(1))*(X(I)-CENTROID(1))
        SUM2=SUM2+(X(I)-CENTROID(1))*(Y(I)-CENTROID(2))
        SUM3=SUM3+(Y(I)-CENTROID(2))*(X(I)-CENTROID(1))
        SUM4=SUM4+(Y(I)-CENTROID(2))*(Y(I)-CENTROID(2))
      ENDDO
      CXX=SUM1/X_POINTS
      CXY=SUM2/X_POINTS
      CYX=SUM3/X_POINTS
      CYY=SUM4/X_POINTS
      COVAR(1,1)=CXX
      COVAR(1,2)=CXY
      COVAR(2,1)=CYX
      COVAR(2,2)=CYY
      LAMBDA1=0.5*(CXX+CYY+SQRT((CXX+CYY)*(CXX+CYY)-4.0*
     1                          (CXX*CYY-CXY*CXY)))
      LAMBDA2=0.5*(CXX+CYY-SQRT((CXX+CYY)*(CXX+CYY)-4.0*
     1                          (CXX*CYY-CXY*CXY)))
      ECCEN=LAMBDA2/LAMBDA1
C
C  Circularity ratio CIRC_RATIO
C
      SUM1=0.0
      DO I=1,N_POINTS
        DI=SQRT((X(I)-CENTROID(1))*(X(I)-CENTROID(1))+
     1          (Y(I)-CENTROID(2))*(Y(I)-CENTROID(2)))
        SUM1=SUM1+DI
      ENDDO
      MUR=SUM1/X_POINTS
      SUM2=0.0
      DO I=1,N_POINTS
        DI=SQRT((X(I)-CENTROID(1))*(X(I)-CENTROID(1))+
     1          (Y(I)-CENTROID(2))*(Y(I)-CENTROID(2)))
        SUM2=SUM2+(DI-MUR)*(DI-MUR)
      ENDDO
      SIGMAR=SQRT(SUM2/X_POINTS)
      CIRC_RATIO=SIGMAR/MUR
C
      IF(SH_AN.EQ.'CDIS') THEN
C
C  Centroid distance function C_DIST(I)
C
        DO I=1,N_POINTS
          C_DIST(I)=SQRT((X(I)-CENTROID(1))*(X(I)-CENTROID(1))+
     1                   (Y(I)-CENTROID(2))*(Y(I)-CENTROID(2)))
        ENDDO
C
      ELSEIF(SH_AN.EQ.'TANG') THEN
C
C  Tangent angle T_ANGLE(I) and normalized tangent angle T2_ANGLE(I)
C
        T_ANGLE(1)=ATAN((Y(1)-Y(N_POINTS))/(X(1)-X(N_POINTS)))*180./PI
        T2_ANGLE(1)=0.
        DO I=2,N_POINTS
          T_ANGLE(I)=ATAN((Y(I)-Y(I-1))/(X(I)-X(I-1)))*180./PI
          T2_ANGLE(I)=T_ANGLE(I)-T_ANGLE(1)-360.0*FLOAT(I)/X_POINTS
        ENDDO
C
      ELSEIF(SH_AN.EQ.'CURV') THEN
C
C  Curvature function K(I) and normalized curvature function K_NORM(I)
C
C......... Computing the first and second derivations of X(I) and Y(I)
C
        STEP=1.
        N_CALC=3
        I_FLAG=2
        CALL DERIV(X,N_POINTS,D1X,D2X,F3,F4,F5,N_CALC,STEP,I_FLAG)
        CALL DERIV(Y,N_POINTS,D1Y,D2Y,F3,F4,F5,N_CALC,STEP,I_FLAG)
C
C......... Computing the curvature function K(I)
C
        SUM1=0.0
        DO I=1,N_POINTS
          NUM=D1X(I)*D2Y(I)-D1Y(I)*D2X(I)
          DEN=(D1X(I)*D1X(I)+D1Y(I)*D1Y(I))**1.5
          K(I)=NUM/DEN
          SUM1=SUM1+ABS(K(I))
        ENDDO
C
C......... Computing the normalized curvature function K_NORM(I)
C
        DO I=1,N_POINTS
          K_NORM(I)=K(I)*X_POINTS/SUM1
        ENDDO
C
C  Average bending energy BE
C
        SUM1=0.0
        DO I=1,N_POINTS
          SUM1=SUM1+K(I)*K(I)
        ENDDO
        BE=SUM1/X_POINTS
C
      ELSEIF(SH_AN.EQ.'TRAR') THEN
C
C  Triangle-area representation TAR(I,J) (normalized to contour area)
C
        NJ=(N_POINTS-1)/2
        DO I=1,N_POINTS
          DO J=1,NJ
            IF((I+J).LE.N_POINTS) THEN
              M=I+J
            ELSE
              M=I+J-N_POINTS
            ENDIF
            IF((I-J).GE.1) THEN
              N=I-J
            ELSE
              N=N_POINTS+I-J
            ENDIF
            TAR(I,J)=(X(N)*Y(I)+Y(N)*X(M)+X(I)*Y(M)-
     1                X(N)*Y(M)-Y(N)*X(I)-Y(I)*X(M))*0.5/C_AREA
          ENDDO
        ENDDO
C
      ELSEIF(SH_AN.EQ.'BEAS') THEN
C
C  Beam angle statistics BAS(I,J)
C
        NJ=(N_POINTS-1)/2
        DO I=1,N_POINTS
          DO J=1,NJ
            IF((I+J).LE.N_POINTS) THEN
              M=I+J
            ELSE
              M=I+J-N_POINTS
            ENDIF
            IF((I-J).GE.1) THEN
              N=I-J
            ELSE
              N=N_POINTS+I-J
            ENDIF
            TH_IPJ=ATAN((Y(M)-Y(I))/(X(M)-X(I)))
            TH_IMJ=ATAN((Y(N)-Y(I))/(X(N)-X(I)))
            BAS(I,J)=(TH_IMJ-TH_IPJ)*180./PI
          ENDDO
        ENDDO
C
      ELSEIF(SH_AN.EQ.'8CCH') THEN
C
C  8-connectivity anticlockwise chain code:
C
C         0: THETA = 0    direction
C         1: THETA = PI/4 direction
C         2: THETA = PI/2
C             ..........
C         7: THETA =-PI/4 direction
C
        DO I=1,N_POINTS
C
C......... Test for THETA = +/- PI/2 (Tangent infinite)
C
          IF(ABS(X(I+1)-X(I)).LT.EPS) THEN
            IF(Y(I+1).GT.Y(I)) THEN
              CHAIN_CODE(I)=2
            ELSE
              CHAIN_CODE(I)=6
            ENDIF
            GOTO 15
          ENDIF
C
          TANGENT=(Y(I+1)-Y(I))/(X(I+1)-X(I))
C
C......... Dividing into 8 half-quadrants
C
          IF((Y(I+1).GE.Y(I)).AND.(X(I+1).GE.X(I))) THEN
C
C.............. Top right quadrant [0,PI/2[
C
            IF(TANGENT.GE.TAN_PIO4) THEN
              IF(TANGENT.GE.TAN_3PIO8) THEN
                CHAIN_CODE(I)=2
              ELSE
                CHAIN_CODE(I)=1
              ENDIF
            ELSE
              IF(TANGENT.GE.TAN_PIO8) THEN
                CHAIN_CODE(I)=1
              ELSE
                CHAIN_CODE(I)=0
              ENDIF
            ENDIF
C
          ELSEIF((Y(I+1).GE.Y(I)).AND.(X(I+1).LT.X(I))) THEN
C
C.............. Top left quadrant ]PI/2,PI]
C
            IF(TANGENT.GE.TAN_3PIO4) THEN
              IF(TANGENT.GE.TAN_7PIO8) THEN
                CHAIN_CODE(I)=4
              ELSE
                CHAIN_CODE(I)=3
              ENDIF
            ELSE
              IF(TANGENT.GE.TAN_5PIO8) THEN
                CHAIN_CODE(I)=3
              ELSE
                CHAIN_CODE(I)=2
              ENDIF
            ENDIF
C
          ELSEIF((Y(I+1).LT.Y(I)).AND.(X(I+1).GE.X(I))) THEN
C
C.............. Bottom right quadrant [0,-PI/2[
C
            IF(TANGENT.GE.TAN_MPIO4) THEN
              IF(TANGENT.GE.TAN_MPIO8) THEN
                CHAIN_CODE(I)=0
              ELSE
                CHAIN_CODE(I)=7
              ENDIF
            ELSE
              IF(TANGENT.GE.TAN_M3PIO8) THEN
                CHAIN_CODE(I)=7
              ELSE
                CHAIN_CODE(I)=6
              ENDIF
            ENDIF
C
          ELSEIF((Y(I+1).LT.Y(I)).AND.(X(I+1).LT.X(I))) THEN
C
C.............. Bottom left quadrant  ]-PI/2,-PI]
C
            IF(TANGENT.GE.TAN_M3PIO4) THEN
              IF(TANGENT.GE.TAN_M5PIO8) THEN
                CHAIN_CODE(I)=6
              ELSE
                CHAIN_CODE(I)=5
              ENDIF
            ELSE
              IF(TANGENT.GE.TAN_M7PIO8) THEN
                CHAIN_CODE(I)=5
              ELSE
                CHAIN_CODE(I)=4
              ENDIF
            ENDIF
          ENDIF
C
  15      CONTINUE
        ENDDO
C
      ELSEIF((SH_AN.EQ.'CLEN').OR.(SH_AN.EQ.'CANG')) THEN
C
C  Chord length and angles distribution
C
        CHORD_DIS_MIN=1.0E+30
        CHORD_DIS_MAX=0.00
        DO I=1,N_POINTS
          CHORD_DIS(I,I)=0.00
          CHORD_ANG(I,I)=0.00
          DO J=1,N_POINTS
C
C......... Computing chord lengths and chord angles
C......... and the min/max
C
            IF(J.EQ.I) GOTO 45
            TMP1=SQRT((X(J)-X(I))*(X(J)-X(I))+
     1                (Y(J)-Y(I))*(Y(J)-Y(I)))
            TMP2=(Y(J)-Y(I))/(X(J)-X(I))
            CHORD_DIS(I,J)=TMP1
            CHORD_DIS_MIN=MIN(TMP1,CHORD_DIS_MIN)
            CHORD_DIS_MAX=MAX(TMP1,CHORD_DIS_MAX)
            CHORD_ANG(I,J)=ATAN(TMP2)*180./PI
  45        CONTINUE
          ENDDO
        ENDDO
C
C......... Computing distributions: putting chords into bins
C
        SIZE_BIN=(CHORD_DIS_MAX-CHORD_DIS_MIN)/FLOAT(N_BIN)
        ANGL_BIN=180./FLOAT(N_BIN)
C
        DO I=1,N_POINTS
          DO J=1,N_POINTS
            IF(J.EQ.I) GOTO 55
C
            IF(SH_AN.EQ.'CLEN') THEN
C
C.............. Chord lengths
C
              CD=CHORD_DIS(I,J)
              DO L=1,N_BIN
                DIS_BIN_LO=CHORD_DIS_MIN+FLOAT(L-1)*SIZE_BIN
                DIS_BIN_UP=DIS_BIN_LO+SIZE_BIN
                JBIN=0
                IF(CD.GE.DIS_BIN_LO) JBIN=JBIN+1
                IF(CD.LE.DIS_BIN_UP) JBIN=JBIN+1
                IF(JBIN.EQ.2) THEN
                  CHORD_DIS_HIS(L)=CHORD_DIS_HIS(L)+1
                  GOTO 25
                ENDIF
              ENDDO
  25          CONTINUE
C
            ELSEIF(SH_AN.EQ.'CANG') THEN
C
C.............. Chord angles
C
              DO L=1,N_BIN
                ANG_BIN_LO=-90.0+FLOAT(L-1)*ANGL_BIN
                ANG_BIN_UP=ANG_BIN_LO+ANGL_BIN
                JBIN=0
                IF(CHORD_ANG(I,J).GE.ANG_BIN_LO) JBIN=JBIN+1
                IF(CHORD_ANG(I,J).LE.ANG_BIN_UP) JBIN=JBIN+1
                IF(JBIN.EQ.2) THEN
                  CHORD_ANG_HIS(L)=CHORD_ANG_HIS(L)+1
                  GOTO 35
                ENDIF
              ENDDO
  35          CONTINUE
C
            ENDIF
C
  55        CONTINUE
C
          ENDDO
        ENDDO
C
      ELSEIF(SH_AN.EQ.'ACDI') THEN
C
          J=N_LEN
          DO I=1,N_POINTS
            IF((I+J).LE.N_POINTS) THEN
              M=I+J
            ELSE
              M=I+J-N_POINTS
            ENDIF
            IF((I-J).GE.1) THEN
              N=I-J
            ELSE
              N=N_POINTS+I-J
            ENDIF
            A=(Y(M)-Y(N))/(X(M)-X(N))
            D_PLUS=ABS(A*(X(M)-X(I))+Y(I)-Y(M))/SQRT(A*A+1.0)
            ACD(I)=D_PLUS
          ENDDO
C
      ELSEIF(SH_AN.EQ.'FOUR') THEN
C
C  Fourier descriptor
C
        DO I=1,N_BIN
          XI=FLOAT(I)
          SUMC=ZEROC
          DO N=1,N_POINTS
            Z(N)=X(N)+IC*Y(N)
            XN=FLOAT(N)
            SUMC=SUMC+Z(N)*CEXP(-TWOPI*IC*XI*XN/X_POINTS)
          ENDDO
          AN(I)=SUMC/SQRT(X_POINTS)
        ENDDO
C
C.............. Computing the inverse transform for checking
C
        DO I=1,N_POINTS
          XI=FLOAT(I)
          SUMC=ZEROC
          DO N=1,N_BIN
            XN=FLOAT(N)
            SUMC=SUMC+AN(N)*CEXP(TWOPI*IC*XI*XN/X_POINTS)
          ENDDO
          XX=REAL(REAL(SUMC))/SQRT(X_POINTS)
          YY=REAL(AIMAG(SUMC))/SQRT(X_POINTS)
          WRITE(97,*) XX,YY
        ENDDO
C
      ENDIF
C
C  Printing general shape information on contour
C
      WRITE(IUO1,20) FILE2
      WRITE(IUO1,30) PERIM
      WRITE(IUO1,40) C_AREA
      WRITE(IUO1,110) COMPACT
      WRITE(IUO1,50) CENTROID(1),CENTROID(2)
      WRITE(IUO1,120) MUR
      WRITE(IUO1,60) THETA
      WRITE(IUO1,61) LAMBDA1
      WRITE(IUO1,62) LAMBDA2
      WRITE(IUO1,70) COVAR(1,1),COVAR(1,2),COVAR(2,1),COVAR(1,2)
      WRITE(IUO1,80) ECCEN
      WRITE(IUO1,90) CIRC_RATIO
      IF(SH_AN.EQ.'CURV') THEN
        WRITE(IUO1,100) BE
      ENDIF
C
C  Writing the shape descriptor into SA
C
      IF(SH_AN.EQ.'CDIS') THEN
        DO I=1,N_POINTS
          SA(I)=C_DIST(I)
          WRITE(98,*) I,SA(I)
        ENDDO
      ELSEIF(SH_AN.EQ.'TANG') THEN
        DO I=1,N_POINTS+1
          SA(I)=T2_ANGLE(I)
          WRITE(98,*) I,SA(I)
        ENDDO
      ELSEIF(SH_AN.EQ.'CURV') THEN
        DO I=1,N_POINTS
          SA(I)=K_NORM(I)
          WRITE(98,*) I,SA(I)
        ENDDO
      ELSEIF(SH_AN.EQ.'TRAR') THEN
        J=N_LEN
        DO I=1,N_POINTS
          SA(I)=TAR(I,J)
          WRITE(98,*) I,SA(I)
        ENDDO
      ELSEIF(SH_AN.EQ.'BEAS') THEN
        J=N_LEN
        DO I=1,N_POINTS
          SA(I)=BAS(I,J)
          WRITE(98,*) I,SA(I)
        ENDDO
      ELSEIF(SH_AN.EQ.'8CCH') THEN
        DO I=1,N_POINTS
          SA(I)=CHAIN_CODE(I)
          WRITE(98,*) I,SA(I)
        ENDDO
      ELSEIF(SH_AN.EQ.'CLEN') THEN 
        IF(N_LEN.EQ.0) THEN
          DO L=1,N_BIN
            SA(L)=CHORD_DIS_HIS(L)
            WRITE(98,*) L,SA(L)
          ENDDO
        ELSE
          J=N_LEN
          DO I=1,N_POINTS
            SA(I)=CHORD_DIS(I,J)
            WRITE(98,*) I,SA(I)
          ENDDO
        ENDIF
      ELSEIF(SH_AN.EQ.'CANG') THEN
        IF(N_LEN.EQ.0) THEN
          DO L=1,N_BIN
            SA(L)=CHORD_ANG_HIS(L)
            WRITE(98,*) L,SA(L)
          ENDDO
        ELSE
          J=N_LEN
          DO I=1,N_POINTS
            SA(I)=CHORD_ANG(I,J)
            WRITE(98,*) I,SA(I)
          ENDDO
        ENDIF
      ELSEIF(SH_AN.EQ.'ACDI') THEN
        DO I=1,N_POINTS
          SA(I)=ACD(I)
          WRITE(98,*) I,SA(I)
        ENDDO
      ELSEIF(SH_AN.EQ.'FOUR') THEN
        DO I=1,N_BIN
          XAI=REAL(REAL(AN(I)))
          YAI=REAL(AIMAG(AN(I)))
          IF(I_FOU.EQ.1) THEN
            SA(I)=SQRT(XAI*XAI+YAI*YAI)
          ELSEIF(I_FOU.EQ.2) THEN
            SA(I)=ATAN2(YAI,XAI)*180./PI
          ELSEIF(I_FOU.EQ.3) THEN
            SA(I)=XAI
          ELSEIF(I_FOU.EQ.4) THEN
            SA(I)=YAI
          ENDIF
          WRITE(98,*) I,SA(I)
        ENDDO
      ENDIF
C
C  Formats
C
  10  FORMAT(//,10X,'<<<<<  N_SIZE SHOULD BE LARGER THAN ',I5,
     1      '  >>>>>',//)
  11  FORMAT(//,10X,'<<<<<  NMAX SHOULD BE LARGER THAN ',I5,
     1      '  >>>>>',//)
  12  FORMAT(//,10X,'<<<<<  NBIN CANNOT BE EQUAL TO ZERO  >>>>>',//)
  13  FORMAT(//,10X,'<<<<<  NBIN TOO LARGE: SET TO N_POINTS  >>>>>',//)
  20  FORMAT(//,6X,'GENERAL SHAPE INFORMATION ON CONTOUR ',A50,/,6X,'|')
  30  FORMAT(6X,'|',3X,'PERIMETER                  : ',E12.6)
  40  FORMAT(6X,'|',3X,'AREA                       : ',E12.6)
  50  FORMAT(6X,'|',3X,'CENTROID                   : ','(',E12.6,',',
     1       E12.6,')')
  60  FORMAT(6X,'|',3X,'AXIS OF LEAST INERTIA ANGLE: ',E12.6)
  61  FORMAT(6X,'|',3X,'LENGTH OF PRINCIPAL AXIS 1 : ',E12.6)
  62  FORMAT(6X,'|',3X,'LENGTH OF PRINCIPAL AXIS 2 : ',E12.6)
  70  FORMAT(6X,'|',3X,'COVARIANCE MATRIX          : (',E12.6,',',
     1       E12.6,')',/,6X,'|',3X,29X,'(',E12.6,',',E12.6,')')
  80  FORMAT(6X,'|',3X,'ECCENTRICITY               : ',E12.6)
  90  FORMAT(6X,'|',3X,'CIRCULARITY RATIO          : ',E12.6)
 100  FORMAT(6X,'|',3X,'AVERAGE BENDING ENERGY     : ',E12.6)
 110  FORMAT(6X,'|',3X,'COMPACTNESS                : ',E12.6)
 120  FORMAT(6X,'|',3X,'MEAN RADIUS                : ',E12.6)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE NORMALIZE_CONTOUR(X,Y,N_POINTS,I_METHOD,IUO1)
C
C  This subroutine normalizes a contour so that different contours
C                        can be compared
C
C  Input parameters:
C
C       * (X,Y)    : points of the input curve
C       * N_POINTS : number of points
C       * I_METHOD : normalization method
C
C                    = 1 : second order moments equal to unity 
C                    = 2 : affine orthogonalization method 
C                    = 3 : equal area method (area = 1) for contour
C                    = 4 : Gu-Kundu method
C
C       * IUO1     : output check file number for printing
C
C
C  Output parameters:
C
C       * (X,Y)    : points of the output curve
C
C  References:  M. Avrithis, Y. Xirouhakis and S. Kolias, 
C               Machine Vision and Applications, 13, 80-94 (2001)
C
C               S. Gu and S. Kundu
C               in proceeding of: Seventh International Conference on 
C               Advances in Pattern Recognition, ICAPR 2009, Kolkata,
C               India, 4-6 February 2009
C
C   Author :  D. Sébilleau
C              
C                                          Last modified : 18 Aug 2014
C
      PARAMETER (N_SIZE=1000)
C
      COMPLEX*8 U1,UN,IC,CSUM1,CSUM2
C
      REAL*4 X(N_SIZE),Y(N_SIZE),X_POINTS
      REAL*4 X_TMP(N_SIZE),Y_TMP(N_SIZE)
      REAL*4 CENTROID(2),C_AREA
      REAL*4 SUM1,SUM2,SUM3,XDIST,YDIST,DIST,MAXD,SIGX,SIGY
      REAL*4 MUX,MUY,SIGMAX,SIGMAY,TAUX,TAUY,M12,M21,RZ
C
      INTEGER P
C
      DATA COEF              /0.707107/
      DATA PI,TWOPI,FOURPI   /3.141593,6.283185,12.566371/
C
      IC=(0.,1.)
C
      IF(N_POINTS.GT.(N_SIZE-1)) THEN
         WRITE(IUO1,10) N_POINTS+1
         STOP
      ENDIF
C
      X_POINTS=FLOAT(N_POINTS)
      X(N_POINTS+1)=X(1)
      Y(N_POINTS+1)=Y(1)
C
      IF(I_METHOD.EQ.1) THEN
C
C.......        Avrithis-Xirouhakis-Kolias method        .......
C.......         at stopped at step 2 (curve S2)         .......
C
C.......  Computation of the moments of the original curve S
C
        SUM1=0.0
        SUM2=0.0
        DO I=1,N_POINTS
          SUM1=SUM1+X(I)
          SUM2=SUM2+Y(I)
        ENDDO
        MUX=SUM1/X_POINTS
        MUY=SUM2/X_POINTS
C
C.......  New curve S1 with center of gravity at origin
C
        DO I=1,N_POINTS
          X_TMP(I)=X(I)-MUX
          Y_TMP(I)=Y(I)-MUY
        ENDDO
C
C.......  Computation of the moments of curve S1
C
        SUM1=0.0
        SUM2=0.0
        DO I=1,N_POINTS
          SUM1=SUM1+X_TMP(I)*X_TMP(I)
          SUM2=SUM2+Y_TMP(I)*Y_TMP(I)
        ENDDO
        SIGMAX=1.0/SQRT(SUM1/X_POINTS)
        SIGMAY=1.0/SQRT(SUM2/X_POINTS)
C
C.......  New curve S2 scaled horizontally and vertically
C
        DO I=1,N_POINTS
          X(I)=X_TMP(I)*SIGMAX
          Y(I)=Y_TMP(I)*SIGMAY
        ENDDO
C
      ELSEIF(I_METHOD.EQ.2) THEN
C
C.......        Avrithis-Xirouhakis-Kolias method        .......
C
C.......  Computation of the moments of the original curve S
C
        SUM1=0.0
        SUM2=0.0
        DO I=1,N_POINTS
          SUM1=SUM1+X(I)
          SUM2=SUM2+Y(I)
        ENDDO
        MUX=SUM1/X_POINTS
        MUY=SUM2/X_POINTS
C
C.......  New curve S1 with center of gravity at origin
C
        DO I=1,N_POINTS
          X_TMP(I)=X(I)-MUX
          Y_TMP(I)=Y(I)-MUY
        ENDDO
C
C.......  Computation of the moments of curve S1
C
        SUM1=0.0
        SUM2=0.0
        DO I=1,N_POINTS
          SUM1=SUM1+X_TMP(I)*X_TMP(I)
          SUM2=SUM2+Y_TMP(I)*Y_TMP(I)
        ENDDO
        SIGMAX=1.0/SQRT(SUM1/X_POINTS)
        SIGMAY=1.0/SQRT(SUM2/X_POINTS)
C
C.......  New curve S2 scaled horizontally and vertically
C
        DO I=1,N_POINTS
          X_TMP(I)=X_TMP(I)*SIGMAX
          Y_TMP(I)=Y_TMP(I)*SIGMAY
        ENDDO
C
C.......  New curve S3 rotated by pi/4
C
        DO I=1,N_POINTS
          X_TMP(I)=COEF*(X_TMP(I)-Y_TMP(I))
          Y_TMP(I)=COEF*(X_TMP(I)+Y_TMP(I))
        ENDDO
C
C.......  Computation of the moments of curve S3
C
        SUM1=0.0
        SUM2=0.0
        DO I=1,N_POINTS
          SUM1=SUM1+X_TMP(I)*X_TMP(I)
          SUM2=SUM2+Y_TMP(I)*Y_TMP(I)
        ENDDO
        TAUX=1.0/SQRT(SUM1/X_POINTS)
        TAUY=1.0/SQRT(SUM2/X_POINTS)
C
C.......  Orthogonalized curve S4
C
        DO I=1,N_POINTS
          X(I)=TAUX*X_TMP(I)
          Y(I)=TAUY*Y_TMP(I)
        ENDDO
C
C.......  Computation of Fourier coefficients U1 and UN  
C.......       for starting point normalization
C
C    As our contour is from 1 to N_POINTS and not 0 to N_POINTS-1,
C         we set X(0)=X(N_POINTS) and Y(0)=Y(N_POINTS) and
C                      go from 0 to N_POINTS-1
C
        CSUM1=(X(N_POINTS)+IC*Y(N_POINTS))
        CSUM2=(X(N_POINTS)+IC*Y(N_POINTS))
        DO I=1,N_POINTS-1
          CSUM1=CSUM1+(X(I)+IC*Y(I))*CEXP(-IC*TWOPI*FLOAT(I)/X_POINTS)
          CSUM2=CSUM2+(X(I)+IC*Y(I))*CEXP(-IC*TWOPI*FLOAT(I)*
     1                (X_POINTS-1)/X_POINTS)
        ENDDO
        U1=CSUM1
        UN=CSUM2
        X1=REAL(U1)
        Y1=REAL(AIMAG(U1))
        XN=REAL(UN)
        YN=REAL(AIMAG(UN))
        A1=ATAN2(Y1,X1)
        AN=ATAN2(YN,XN)
        P=INT(X_POINTS*(A1-AN)/FOURPI)
        P=MOD(P,N_POINTS/2)
C
C.......  Shifting the contour by -P for starting point normalization
C
        DO I=1,N_POINTS
          IF((I-P).LT.1) THEN
            X_TMP(I)=X(N_POINTS+I-P)
            Y_TMP(I)=Y(N_POINTS+I-P)
          ELSE
            X_TMP(I)=X(I-P)
            Y_TMP(I)=Y(I-P)
          ENDIF
        ENDDO
        DO I=1,N_POINTS
          X(I)=X_TMP(I)
          Y(I)=Y_TMP(I)
        ENDDO
C
C.......  Rotation and reflection normalization  (curve Z2)
C
C
C.......  Computation of Fourier coefficients U1 and UN  
C.......      of starting point normalized contour (contour Z)
C
        CSUM1=(X(N_POINTS)+IC*Y(N_POINTS))
        CSUM2=(X(N_POINTS)+IC*Y(N_POINTS))
        DO I=1,N_POINTS-1
          CSUM1=CSUM1+(X(I)+IC*Y(I))*CEXP(-IC*TWOPI*FLOAT(I)/X_POINTS)
          CSUM2=CSUM2+(X(I)+IC*Y(I))*CEXP(-IC*TWOPI*FLOAT(I)*
     1                (X_POINTS-1)/X_POINTS)
        ENDDO
        U1=CSUM1
        UN=CSUM2
        X1=REAL(U1)
        Y1=REAL(AIMAG(U1))
        XN=REAL(UN)
        YN=REAL(AIMAG(UN))
        A1=ATAN2(Y1,X1)
        AN=ATAN2(YN,XN)
        RZ=MOD(0.5*(A1+AN),PI)
C
C.......  Computation of new contour Z1  and its moments M12 and M21
C
        DO I=1,N_POINTS
          X(I)=X(I)*COS(RZ)-Y(I)*SIN(RZ)
          Y(I)=Y(I)*COS(RZ)+X(I)*SIN(RZ)
        ENDDO
C
        SUM1=0.0
        SUM2=0.0
        DO I=1,N_POINTS
          SUM1=SUM1+X(I)*Y(I)*Y(I)
          SUM2=SUM2+X(I)*X(I)*Y(I)
        ENDDO
        M12=SUM1/X_POINTS
        M21=SUM2/X_POINTS
C
        IF(M12.GE.0.0) THEN
          SIGX=1.0
        ELSE
          SIGX=-1.0
        ENDIF
        IF(M21.GE.0.0) THEN
          SIGY=1.0
        ELSE
          SIGY=-1.0
        ENDIF
C
C.......  Contour Z2
C
        DO I=1,N_POINTS
          X(I)=SIGX*X(I)
          Y(I)=SIGY*Y(I)
        ENDDO
C
      ELSEIF(I_METHOD.EQ.3) THEN
C
C.......        Simple normalization to unit area        .......
C
C.......  Computation of the centroid and area of contour
C
        SUM1=0.0
        DO I=1,N_POINTS
           SUM1=SUM1+(X(I)*Y(I+1)-X(I+1)*Y(I))
        ENDDO
        C_AREA=0.5*ABS(SUM1)
C
        SUM1=0.0
        SUM2=0.0
        DO I=1,N_POINTS
          SUM1=SUM1+(X(I)+X(I+1))*(X(I)*Y(I+1)-X(I+1)*Y(I))
          SUM2=SUM2+(Y(I)+Y(I+1))*(X(I)*Y(I+1)-X(I+1)*Y(I))
        ENDDO
        CENTROID(1)=SUM1/(6.*C_AREA)
        CENTROID(2)=SUM2/(6.*C_AREA)
C
C.......  Scaling (X,Y) to area and centering on centroid
C
        DO I=1,N_POINTS
          X(I)=(X(I)-CENTROID(1))/SQRT(C_AREA)
          Y(I)=(Y(I)-CENTROID(2))/SQRT(C_AREA)
        ENDDO
C
      ELSEIF(I_METHOD.EQ.4) THEN
C
C.......        Gu-Kundu method        .......
C
C.......  Computation of the centroid 
C
        SUM1=0.0
        SUM2=0.0
        SUM3=0.0
        DO I=1,N_POINTS
          XDIST=(X(I+1)-X(I))*(X(I+1)-X(I))
          YDIST=(Y(I+1)-Y(I))*(Y(I+1)-Y(I))
          DIST=SQRT(XDIST+YDIST)
          SUM1=SUM1+0.5*(X(I+1)+X(I))*DIST
          SUM2=SUM2+0.5*(Y(I+1)+Y(I))*DIST
          SUM3=SUM3+DIST
        ENDDO
        CENTROID(1)=SUM1/SUM3
        CENTROID(2)=SUM2/SUM3
C
C.......  Maximum distance to centroid 
C
        MAXD=0.0
        DO I=1,N_POINTS
          XDIST=(X(I)-CENTROID(1))*(X(I)-CENTROID(1))
          YDIST=(Y(I)-CENTROID(2))*(Y(I)-CENTROID(2))
          DIST=SQRT(XDIST+YDIST)
          MAXD=MAX(DIST,MAXD)
        ENDDO
C
C.......  New contour
C       
        DO I=1,N_POINTS
          X(I)=(X(I)-CENTROID(1))/MAXD
          Y(I)=(Y(I)-CENTROID(2))/MAXD
        ENDDO
C
      ENDIF
C
C  Formats
C
  10  FORMAT(//,10X,'<<<<< DIMENSION ERROR: N_SIZE SHOULD BE >>>>>',/,
     1          10X,'<<<<<   ',I5,' IN ROUTINE NORMALIZE    >>>>>',//)
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE NORMALIZE_CURVE(X,Y,N_POINTS,I_METHOD,IUO1)
C
C  This subroutine normalizes a curve Y=f(X) so that different curves
C                       can be compared
C
C  Input parameters:
C
C       * (X,Y)    : points of the input curve
C       * N_POINTS : number of points
C       * I_METHOD : normalization method
C
C                    = 1 : second order central moment equals to unity 
C                    = 2 : equal area method (area = 1) for curve
C                    = 3 : normalization to maximum
C                    = 4 : decimal scaling
C                    = 5 : normalization/rescaling in [0,1]
C
C       * IUO1     : output check file number for printing
C
C
C  Output parameters:
C
C       * (X,Y)    : points of the output curve
C
C   Author :  D. Sébilleau
C
C              
C                                          Last modified : 18 Aug 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 X(N_SIZE),Y(N_SIZE),X_POINTS
      REAL*4 Y_TMP(N_SIZE)
      REAL*4 SUM1,MUY,SIGMAY
      REAL*4 MAX_PEAK,MIN_PEAK,DIFF_PEAK,DEC_SCALE,C_AREA
      REAL*4 MAXY,TMP,COMP
C
      DATA EPS /0.001/
C
      IF(N_POINTS.GT.(N_SIZE-1)) THEN
         WRITE(IUO1,10) N_POINTS+1
         STOP
      ENDIF
C
      X_POINTS=FLOAT(N_POINTS)
      X(N_POINTS+1)=X(1)
      Y(N_POINTS+1)=Y(1)
C
      IF(I_METHOD.EQ.1) THEN
C
C.......  Normalization so that central second order moment is equal to 1
C.......  (also called z-score normalization)
C
C
C.......  Computation of the moments of the original curve S
C
        SUM1=0.0
        DO I=1,N_POINTS
          SUM1=SUM1+Y(I)
        ENDDO
        MUY=SUM1/X_POINTS
C
C.......  New curve S1 with center of gravity at origin
C
        DO I=1,N_POINTS
          Y_TMP(I)=Y(I)-MUY
        ENDDO
C
C.......  Computation of the moments of curve S1
C
        SUM1=0.0
        DO I=1,N_POINTS
          SUM1=SUM1+Y_TMP(I)*Y_TMP(I)
        ENDDO
c        SIGMAY=SQRT(SUM1/X_POINTS)
        SIGMAY=SQRT(SUM1)
C
C.......  New curve S2 scaled horizontally and vertically
C
        DO I=1,N_POINTS
          Y(I)=Y_TMP(I)/SIGMAY
        ENDDO
C
      ELSEIF(I_METHOD.EQ.2) THEN
C
C.......        Normalization to unit area        .......
C
        SUM1=0.0
        DO I=1,N_POINTS
           SUM1=SUM1+0.5*(Y(I+1)-Y(I))*(X(I+1)-X(I))
        ENDDO
        C_AREA=ABS(SUM1)
C
C.......  Scaling (X,Y) to area 
C
        DO I=1,N_POINTS
          Y(I)=Y(I)/SQRT(C_AREA)
        ENDDO
C
      ELSEIF(I_METHOD.EQ.3) THEN
C
C.......        Normalization to highest peak        .......
C
        MAX_PEAK=-1.0E+8
        DO I=1,N_POINTS
          MAX_PEAK=MAX(MAX_PEAK,Y(I))
        ENDDO
C
        IF(ABS(MAX_PEAK+1.0E+8).LT.EPS) THEN
         WRITE(IUO1,20) 
         STOP
        ENDIF
C
C.......  Scaling (X,Y) to highest peak 
C
        DO I=1,N_POINTS
          Y(I)=Y(I)/MAX_PEAK
        ENDDO
C
      ELSEIF(I_METHOD.EQ.4) THEN
C
C.......        Decimal scaling        .......
C
        MAX_PEAK=-1.0E+8
        DO I=1,N_POINTS
          MAX_PEAK=MAX(MAX_PEAK,Y(I))
        ENDDO
C
        IF(ABS(MAX_PEAK+1.0E+8).LT.EPS) THEN
         WRITE(IUO1,20) 
         STOP
        ENDIF
C
        DO K=0,100
          COMP=10.0**K
          MAXY=-1.0E+8
          DO I=1,N_POINTS
            TMP=ABS(Y(I)/COMP)
            MAXY=MAX(MAXY,TMP)
          ENDDO
          IF(MAXY.LT.1.0) GOTO 15
        ENDDO
C
  15    DEC_SCALE=10.0**K
C
        DO I=1,N_POINTS
          Y(I)=Y(I)/DEC_SCALE
        ENDDO
C
      ELSEIF(I_METHOD.EQ.5) THEN
C
C.......        Normalization to the range [0,1]        .......
C.......        (also called min-max normalization)     .......
C
        MAX_PEAK=-1.0E+8
        MIN_PEAK=1.0E+8
        DO I=1,N_POINTS
          MAX_PEAK=MAX(MAX_PEAK,Y(I))
          MIN_PEAK=MIN(MIN_PEAK,Y(I))
        ENDDO
C
        IF(ABS(MAX_PEAK+1.0E+8).LT.EPS) THEN
         WRITE(IUO1,20) 
         STOP
        ENDIF
        IF(ABS(MIN_PEAK-1.0E+8).LT.EPS) THEN
         WRITE(IUO1,30) 
         STOP
        ENDIF
C
C.......  Scaling (X,Y) in [0,1] 
C
        DIFF_PEAK=MAX_PEAK-MIN_PEAK
        DO I=1,N_POINTS
          Y(I)=(Y(I)-MIN_PEAK)/DIFF_PEAK
        ENDDO
C
      ENDIF
C
C
C  Formats
C
  10  FORMAT(//,10X,'<<<<< DIMENSION ERROR: N_SIZE SHOULD BE >>>>>',/,
     1          10X,'<<<<<   ',I5,' IN ROUTINE NORMALIZE    >>>>>',//)
  20  FORMAT(//,10X,'<<<<<  ERROR IN ROUTINE NORMALIZE_CURVE >>>>>',/,
     1          10X,'<<<<<   MAXPEAK SHOULD BE INITIALIZED   >>>>>',/,
     2          10X,'<<<<<   TO A VALUE LOWER THAN -10E+08   >>>>>'//)
  30  FORMAT(//,10X,'<<<<<  ERROR IN ROUTINE NORMALIZE_CURVE >>>>>',/,
     1          10X,'<<<<<   MINPEAK SHOULD BE INITIALIZED   >>>>>',/,
     2          10X,'<<<<<   TO A VALUE HIGHER THAN 10E+08   >>>>>'//)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE NORMALIZE_COEF(X,EXPE,CALC,N_POINTS,NFILE,I_NORM,IUO1,
     1                          CALCFILE,CNORM)
C
C  This subroutine computes the normalization coefficient 
C    to scale the calculated curve to the experimental curve.
C
C  Both curves are supposed to have the same range and number
C    of points
C
C  Input parameters:
C
C       * EXPE     : points of the experimental curve
C       * CALC     : points of the calculated curve    (I_NORM > 0)
C       * N_POINTS : number of points
C       * NFILE    : number of calculation files       (I_NORM < 0)
C       * I_NORM   : normalization method
C
C                    =  0 : no normalization
C                    =  1 : C = sum |EXPE(I)| / sum |CALC(I)|
C                    =  2 : C = sum |EXPE(I)*CALC(I)| / sum |CALC(I)^2|
C                    = -1 : same as 1 but average of all CALC used
C                    = -2 : same as 2 but average of all CALC used 
C
C       * IUO1     : output check file number for printing
C       * CALCFILE : name of the calculated files      (I_NORM < 0)
C
C  Output parameter:
C
C       * CNORM    : normalization coefficient 
C
C   Author :  D. Sébilleau
C              
C                                          Last modified : 19 Aug 2014
C
      PARAMETER (N_SIZE=1000,N_FILES=1000)
C
      REAL*4 EXPE(N_SIZE),CALC(N_SIZE),Y_MEAN(N_SIZE)
      real*4 x(N_SIZE)
      REAL*4 SUM1,SUM2,CNORM,X_FILES,DUMMY,Y
C
      CHARACTER*40 CALCFILE(N_FILES),FIRSTFILE,AVERFILE
C
      IF(I_NORM.GT.0) THEN
C
C  Each calculation is normalized separately to experiment
C
        SUM1=0.
        SUM2=0.
C
        DO J=1,N_POINTS
C
          IF(I_NORM.EQ.1) THEN
            SUM1=SUM1+ABS(EXPE(J))
            SUM2=SUM2+ABS(CALC(J))
          ELSEIF(I_NORM.EQ.2) THEN
            SUM1=SUM1+ABS(EXPE(J)*CALC(J))
            SUM2=SUM2+CALC(J)*CALC(J)
          ENDIF
C
        ENDDO
C
      ELSEIF(I_NORM.EQ.0) THEN
C
C  No normalization
C
        SUM1=1.0
        SUM2=1.0
C
      ELSEIF(I_NORM.LT.0) THEN
C
C  The average of the calculations normalized to experiment
C
        IF(NFILE.GT.N_FILES) THEN
          WRITE(IUO1,11) NFILE
          STOP
        ENDIF
C
C........  Checking for the name of the directory 
C........  where the calculation files are stored
C
        FIRSTFILE=CALCFILE(1)
        N_SL=1
        DO J_CHAR=1,40
          IF(FIRSTFILE(J_CHAR:J_CHAR).EQ.'/') GOTO 30
          N_SL=N_SL+1
        ENDDO
   30   CONTINUE
        AVERFILE=FIRSTFILE(1:N_SL-1)//'/calculation_ave.dat'
        OPEN(UNIT=55, FILE=AVERFILE, STATUS='unknown')
C
        X_FILES=FLOAT(NFILE)
C
        DO J=1,N_POINTS
          Y_MEAN(J)=0.
        ENDDO
C
C........  Computing the average curve Y_MEAN of the NFILE calculation files
C
        DO JFILE=1,NFILE
C
          NUNIT3=50
          OPEN(UNIT=NUNIT3, FILE=CALCFILE(JFILE), STATUS='unknown')
          DO JLINE=1,N_POINTS
            READ(NUNIT3,*) DUMMY,Y
            Y_MEAN(JLINE)=Y_MEAN(JLINE)+Y/X_FILES
          ENDDO
          CLOSE(NUNIT3)
C
        ENDDO
C
C........  Storage of the averaged calculation file
C
        DO I=1,N_POINTS
          WRITE(55,*) X(I),Y_MEAN(I)
        ENDDO
        CLOSE(55)
C
        SUM1=0.
        SUM2=0.
C
        DO J=1,N_POINTS
C
          IF(I_NORM.EQ.-1) THEN
            SUM1=SUM1+ABS(EXPE(J))
            SUM2=SUM2+ABS(Y_MEAN(J))
          ELSEIF(I_NORM.EQ.-2) THEN
            SUM1=SUM1+ABS(EXPE(J)*Y_MEAN(J))
            SUM2=SUM2+Y_MEAN(J)*Y_MEAN(J)
          ENDIF
C
        ENDDO
C
      ENDIF
C
      CNORM=SUM1/SUM2
C
C  Format
C
  11  FORMAT(//,10X,'<<<<<  N_FILES SHOULD BE LARGER THAN ',I5,
     1      '  >>>>>',//)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE RESCALE_TO_EXP(CALC,N_POINTS,EXP_MIN,EXP_MAX)
C
C  This subroutine rescales the calculated file to the min and max 
C  of the experimental file around the Y = (EXP_MIN+EXP_MAX)/2 value
C
C  Input parameters:
C
C       * CALC     : points of the calculation curve
C       * N_POINTS : number of points in the curves
C       * EXP_MIN  : minimum of the experimental curve
C       * EXP_MAX  : maximum of the experimental curve
C
C  Output parameters:
C
C       * CALC     : points of the rescaled calculation curve
C
C   Author :  D. Sébilleau
C  
C                                           Last modified : 20 Aug 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 CALC(N_SIZE)
      REAL*4 CAL_MIN,CAL_MAX,MEAN_CAL,TMP_CAL
      REAL*4 EXP_MIN,EXP_MAX,MEAN_EXP
      REAL*4 LARGE
C
      DATA LARGE /1.E+30/
C
C
C  Calculation of the maximum, minimum  and mean of calculation
C            within these intersection bounds
C
      CAL_MIN=LARGE
      CAL_MAX=0.
      MEAN_CAL=0.
C
      MEAN_EXP=0.5*(EXP_MIN+EXP_MAX)
C
      DO J=1,N_POINTS
        CAL_MIN=MIN(CAL_MIN,CALC(J))
        CAL_MAX=MAX(CAL_MAX,CALC(J))
      ENDDO
C
      MEAN_CAL=0.5*(CAL_MAX+CAL_MIN)
      CAL_MAX=CAL_MAX-MEAN_CAL
      CAL_MIN=MEAN_CAL-CAL_MIN
C
C  Shifting the calculation to mean = 0 and scaling positive 
C  by (EXP_MAX-MEAN_EXP)/CAL_MAX and negative values  by 
C  (EXP_MIN-MEAN_EXP)/CAL_MIN, and then shifting it back 
C  to MEAN_EXP
C
      DO J=1,N_POINTS
C
        TMP_CAL=CALC(J)-MEAN_CAL
        IF(TMP_CAL.GE.0.) THEN
          TMP_CAL=TMP_CAL*ABS(EXP_MAX-MEAN_EXP)/CAL_MAX
        ELSE
          TMP_CAL=TMP_CAL*ABS(EXP_MIN-MEAN_EXP)/CAL_MIN
        ENDIF
C
        CALC(J)=TMP_CAL+MEAN_EXP
C
      ENDDO
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE SHIFT_CURVE(Y,N_POINTS,SHIFT)
C
C  This subroutine shifts a curve by SHIFT so that none of their values
C     is negative. This is necessary for some comparison methods
C
C  Input parameters:
C
C       * Y        : points of the original curve
C       * N_POINTS : number of points in the curves
C       * SHIFT    : shift value
C
C  Output parameter:
C
C       * Y        : points of the shifted curve
C
C   Author :  D. Sébilleau
C  
C                                           Last modified : 20 Aug 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 Y(N_SIZE)
      REAL*4 SHIFT
C
      DO I=1,N_POINTS
C
        Y(I)=Y(I)+SHIFT
C
      ENDDO
C
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE WEIGHTS(X,EXPE,CALC,N_POINTS,CNORM,IUO1,W)
C
C  This subroutine computes weights for the comparison of two curves
C
C  Input parameters:
C
C       * X        : Abscissa of the curve
C       * EXPE     : points of the experimental curve
C       * CALC     : points of calculated curve
C       * N_POINTS : number of points in the curves
C       * CNORM    : normalization coefficient for calculated curve
C       * IUO1     : output check file number for printing
C
C       * I_WEIGHT : switch to select the weights
C       * ALPHA    :  |
C       * BETA     :  |  parameters for the weights
C       * SIGMA    :  |
C       * I_SHIFT  : switch to introduce a shift in I (I_WEIGHT = 6 and 7)
C       * MAXW     : value of X(I) for which W(I) is maximal
C
C  Output parameter:
C
C       * W        : normalized weights for the curve points
C
C  References:
C               R. De Gelder, R. Wehrens and J. A Hageman,
C               J. Comput. Chem. 22, 273-289 (2001) (I_WEIGHT = 6 and 7)
C
C   Author :  D. Sébilleau
C  
C                                           Last modified : 28 Aug 2014
C
      PARAMETER (N_SIZE=1000)
C
      REAL*4 X(N_SIZE),EXPE(N_SIZE),CALC(N_SIZE),W(N_SIZE)
      REAL*4 SIGMAE(N_POINTS),COSINE(N_POINTS)
      REAL*4 XJ,X_POINTS,ALPHA,BETA,SIGMA,MAXW,MEANE,SUMW,MAXEC,TMP
C
      COMMON /PAR_WEI/   I_WEIGHT,I_SHIFT,ALPHA,BETA,SIGMA,MAXW 
C
      DATA PI,EPS /3.141593,1.0E-5/
C
      X_POINTS=FLOAT(N_POINTS)
C
C  Preliminary calculations for some weights
C
      IF(I_WEIGHT.EQ.3) THEN
C
C..............  Mean value  
C 
        MEANE=0.
        DO J=1,N_POINTS
          MEANE=MEANE+EXPE(J)
        ENDDO
        MEANE=MEANE/X_POINTS
C
C..............  Sigma square
C
        DO J=1,N_POINTS
          SIGMAE(J)=(EXPE(J)-MEANE)*(EXPE(J)-MEANE)
        ENDDO
      ENDIF
C
C  Calculation of the shift index JM (I_WEIGHT = 6, 7 and 8)
C    and of X = 0. index J0
C
      IF(I_SHIFT.GT.0) THEN
        CALL LOCATE(X,N_POINTS,MAXW,JM,1)
        CALL LOCATE(X,N_POINTS,0.,J0,1)
      ELSE
        JM=0
        J0=0
      ENDIF
C
      IF(I_WEIGHT.EQ.5) THEN
        MAXEC=0.00
        DO J=1,N_POINTS
          TMP=ABS(EXPE(J)-CNORM*CALC(J))
          MAXEC=MAX(MAXEC,TMP)
        ENDDO
      ENDIF
C
      IF(I_WEIGHT.EQ.8) THEN
        DO J=1,N_POINTS
          K=J-JM+J0
          IF(K.LE.0) K=J0+JM-J
          IF(K.GT.0) THEN
            COSINE(J)=COS(X(K)*PI/180.)
          ELSE
            COSINE(J)=1.0
          ENDIF
        ENDDO
      ENDIF
C
C  Computing the weights
C
      SUMW=0.
      DO J=1,N_POINTS
C
        IF(I_WEIGHT.EQ.0) THEN
          W(J)=X_POINTS
        ELSEIF(I_WEIGHT.EQ.1) THEN
          W(J)=(EXPE(J)+CNORM*CALC(J))/(2.0*CNORM*EXPE(J)*CALC(J))
        ELSEIF(I_WEIGHT.EQ.2) THEN
          W(J)=(EXPE(J)*EXPE(J)+CNORM*CNORM*CALC(J)*CALC(J))/
     1         (2.0*CNORM*CNORM*EXPE(J)*EXPE(J)*CALC(J)*CALC(J))
        ELSEIF(I_WEIGHT.EQ.3) THEN
          W(J)=1.0/SIGMAE(J)
        ELSEIF(I_WEIGHT.EQ.4) THEN
          W(J)=EXP(-(EXPE(J)-CNORM*CALC(J))*(EXPE(J)-CNORM*CALC(J))/
     1             2.0*SIGMA*SIGMA)/(2.0*PI*SIGMA*SIGMA)
        ELSEIF(I_WEIGHT.EQ.5) THEN
          W(J)=1.0-ABS(EXPE(J)-CNORM*CALC(J))/MAXEC
        ELSEIF(I_WEIGHT.EQ.6) THEN
          W(J)=1.0/(1.0+ALPHA*(ABS(FLOAT(J-JM))**BETA))
        ELSEIF(I_WEIGHT.EQ.7) THEN
          XJ=ABS(FLOAT(J-JM))
          IF(XJ.LT.BETA) THEN
            W(J)=1.0-ALPHA*XJ/BETA
            IF(W(J).LT.EPS) W(J)=0.0
          ELSE
            W(J)=0.0
          ENDIF
        ELSEIF(I_WEIGHT.EQ.8) THEN
          W(J)=EXP(-ALPHA/(COSINE(J)**BETA))
        ENDIF
        SUMW=SUMW+W(J)
C
      ENDDO
C
C  Normalizing the weights to N_POINTS
C
      IF(ABS(SUMW).GT.EPS) THEN
        DO J=1,N_POINTS
          W(J)=W(J)*X_POINTS/SUMW
        ENDDO
      ELSE
        DO J=1,N_POINTS
          W(J)=1.0
        ENDDO
        WRITE(IUO1,10)
      ENDIF
C
C  Format
C
  10  FORMAT(//,10X,'<<<<<  ERROR IN THE CALCULATION OF ',
     1             'THE WEIGHTS:  >>>>>',/,
     2          10X'<<<<<      WJ SET TO 1.0 IN THE CALCULATION',
     3      '      >>>>>',//)
C
      RETURN
C
      END
C
C=======================================================================
C
C                         LAPACK inversion subroutines
C
C=======================================================================
C
C
C======================================================================
C
      SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETRI computes the inverse of a matrix using the LU factorization
*  computed by ZGETRF.
*
*  This method inverts U and then computes inv(A) by solving the system
*  inv(A)*L = inv(U) for inv(A).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the factors L and U from the factorization
*          A = P*L*U as computed by ZGETRF.
*          On exit, if INFO = 0, the inverse of the original matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimal performance LWORK >= N*NB, where NB is
*          the optimal blocksize returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*                singular and its inverse could not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB,
     $                   NBMIN, NN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGEMV, ZSWAP, ZTRSM, ZTRTRI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NB = ILAENV( 1, 'ZGETRI', ' ', N, -1, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRI', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Form inv(U).  If INFO > 0 from ZTRTRI, then U is singular,
*     and the inverse is not computed.
*
      CALL ZTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
      IF( INFO.GT.0 )
     $   RETURN
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = MAX( LDWORK*NB, 1 )
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZGETRI', ' ', N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = N
      END IF
*
*     Solve the equation inv(A)*L = inv(U) for inv(A).
*
      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         DO 20 J = N, 1, -1
*
*           Copy current column of L to WORK and replace with zeros.
*
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
*
*           Compute current column of inv(A).
*
            IF( J.LT.N )
     $         CALL ZGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE
*
*        Use blocked code.
*
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
*
*           Copy current block column of L to WORK and replace with
*           zeros.
*
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
*
*           Compute current block column of inv(A).
*
            IF( J+JB.LE.N )
     $         CALL ZGEMM( 'No transpose', 'No transpose', N, JB,
     $                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
     $                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL ZTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
   50    CONTINUE
      END IF
*
*     Apply column interchanges.
*
      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         IF( JP.NE.J )
     $      CALL ZSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of ZGETRI
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZTRTI2 computes the inverse of a complex upper or lower triangular
*  matrix.
*
*  This is the Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading n by n upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      COMPLEX*16         AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZSCAL, ZTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTRTI2', -INFO )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Compute inverse of upper triangular matrix.
*
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
*
*           Compute elements 1:j-1 of j-th column.
*
            CALL ZTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,
     $                  A( 1, J ), 1 )
            CALL ZSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        Compute inverse of lower triangular matrix.
*
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
*
*              Compute elements j+1:n of j-th column.
*
               CALL ZTRMV( 'Lower', 'No transpose', DIAG, N-J,
     $                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL ZSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of ZTRTI2
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZTRTRI( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZTRTRI computes the inverse of a complex upper or lower triangular
*  matrix A.
*
*  This is the Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*               matrix is singular and its inverse can not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J, JB, NB, NN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZTRMM, ZTRSM, ZTRTI2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTRTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity if non-unit.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
         INFO = 0
      END IF
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'ZTRTRI', UPLO // DIAG, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code
*
         IF( UPPER ) THEN
*
*           Compute inverse of upper triangular matrix
*
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
*
*              Compute rows 1:j-1 of current block column
*
               CALL ZTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, ONE, A, LDA, A( 1, J ), LDA )
               CALL ZTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
*
*              Compute inverse of current diagonal block
*
               CALL ZTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
   20       CONTINUE
         ELSE
*
*           Compute inverse of lower triangular matrix
*
            NN = ( ( N-1 ) / NB )*NB + 1
            DO 30 J = NN, 1, -NB
               JB = MIN( NB, N-J+1 )
               IF( J+JB.LE.N ) THEN
*
*                 Compute rows j+jb:n of current block column
*
                  CALL ZTRMM( 'Left', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
     $                        A( J+JB, J ), LDA )
                  CALL ZTRSM( 'Right', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, -ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
*
*              Compute inverse of current diagonal block
*
               CALL ZTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
   30       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZTRTRI
*
      END
C
C
C======================================================================
C
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  Purpose
*  =======
*
*  IEEECK is called from the ILAENV to verify that Infinity and
*  possibly NaN arithmetic is safe (i.e. will not trap).
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies whether to test just for inifinity arithmetic
*          or whether to test for infinity and NaN arithmetic.
*          = 0: Verify infinity arithmetic only.
*          = 1: Verify infinity and NaN arithmetic.
*
*  ZERO    (input) REAL
*          Must contain the value 0.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  ONE     (input) REAL
*          Must contain the value 1.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  RETURN VALUE:  INTEGER
*          = 0:  Arithmetic failed to produce the correct answers
*          = 1:  Arithmetic produced the correct answers
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. Executable Statements ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      RETURN
      END
C
C
C======================================================================
C
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*  -- LAPACK auxiliary routine (version 3.1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     January 2007
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  ILAENV returns an INTEGER
*  if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
*  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines (DEPRECATED)
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR method
*               for nonsymmetric eigenvalue problems (DEPRECATED)
*          = 9: maximum size of the subproblems at the bottom of the
*               computation tree in the divide-and-conquer algorithm
*               (used by xGELSD and xGESDD)
*          =10: ieee NaN arithmetic can be trusted not to trap
*          =11: infinity arithmetic can be trusted not to trap
*          12 <= ISPEC <= 16:
*               xHSEQR or one of its subroutines,
*               see IPARMQ for detailed explanation
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER            IEEECK, IPARMQ
      EXTERNAL           IEEECK, IPARMQ
*     ..
*     .. Executable Statements ..
*
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,
     $        130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
   10 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $             I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
*
      GO TO ( 50, 60, 70 )ISPEC
*
   50 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
   60 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
   70 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
   80 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
   90 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  100 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  110 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  120 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  130 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
  140 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
*
  150 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
*
  160 CONTINUE
*
*     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. 
*
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
*
*     End of ILAENV
*
      END
C
C
C======================================================================
C
      LOGICAL FUNCTION LSAME(CA,CB)
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER CA,CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER INTA,INTB,ZCODE
*     ..
*
*     Test if the characters are equal
*
      LSAME = CA .EQ. CB
      IF (LSAME) RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR('Z')
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
*
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
          IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
          IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
*
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
          IF (INTA.GE.129 .AND. INTA.LE.137 .OR.
     +        INTA.GE.145 .AND. INTA.LE.153 .OR.
     +        INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
          IF (INTB.GE.129 .AND. INTB.LE.137 .OR.
     +        INTB.GE.145 .AND. INTB.LE.153 .OR.
     +        INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
*
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
          IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
          IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
      LSAME = INTA .EQ. INTB
*
*     RETURN
*
*     End of LSAME
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN
      INTEGER            I, J, JP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            IZAMAX
      EXTERNAL           DLAMCH, IZAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Compute machine safe minimum
*
      SFMIN = DLAMCH('S') 
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IZAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL ZSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M ) THEN
               IF( ABS(A( J, J )) .GE. SFMIN ) THEN
                  CALL ZSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
               ELSE
                  DO 20 I = 1, M-J
                     A( J+I, J ) = A( J+I, J ) / A( J, J )
   20             CONTINUE
               END IF
            END IF
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL ZGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),
     $                  LDA, A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of ZGETF2
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGETF2, ZLASWP, ZTRSM
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'ZGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL ZGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL ZGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL ZLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL ZLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL ZGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZGETRF
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
*  Further Details
*  ===============
*
*  Modified by
*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      COMPLEX*16         TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASWP
*
      END
C
C
C======================================================================
C
      SUBROUTINE XERBLA(SRNAME,INFO)
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER INFO
      CHARACTER*6 SRNAME
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*
      WRITE (*,FMT=9999) SRNAME,INFO
*
      STOP
*
 9999 FORMAT (' ** On entry to ',A6,' parameter number ',I2,' had ',
     +       'an illegal value')
*
*     End of XERBLA
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  Purpose
*  =======
*
*  ZGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
*     ..
*     .. Local Scalars ..
      DOUBLE COMPLEX TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL CONJA,CONJB,NOTA,NOTB
*     ..
*     .. Parameters ..
      DOUBLE COMPLEX ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      DOUBLE COMPLEX ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
*     B  respectively are to be  transposed but  not conjugated  and set
*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
*     and the number of rows of  B  respectively.
*
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      CONJA = LSAME(TRANSA,'C')
      CONJB = LSAME(TRANSB,'C')
      IF (NOTA) THEN
          NROWA = M
          NCOLA = K
      ELSE
          NROWA = K
          NCOLA = M
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.CONJA) .AND.
     +    (.NOT.LSAME(TRANSA,'T'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.CONJB) .AND.
     +         (.NOT.LSAME(TRANSB,'T'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGEMM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (NOTB) THEN
          IF (NOTA) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
              DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                      IF (B(L,J).NE.ZERO) THEN
                          TEMP = ALPHA*B(L,J)
                          DO 70 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
   70                     CONTINUE
                      END IF
   80             CONTINUE
   90         CONTINUE
          ELSE IF (CONJA) THEN
*
*           Form  C := alpha*conjg( A' )*B + beta*C.
*
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120         CONTINUE
          ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
              DO 150 J = 1,N
                  DO 140 I = 1,M
                      TEMP = ZERO
                      DO 130 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  130                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  140             CONTINUE
  150         CONTINUE
          END IF
      ELSE IF (NOTA) THEN
          IF (CONJB) THEN
*
*           Form  C := alpha*A*conjg( B' ) + beta*C.
*
              DO 200 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 160 I = 1,M
                          C(I,J) = ZERO
  160                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 170 I = 1,M
                          C(I,J) = BETA*C(I,J)
  170                 CONTINUE
                  END IF
                  DO 190 L = 1,K
                      IF (B(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*DCONJG(B(J,L))
                          DO 180 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  180                     CONTINUE
                      END IF
  190             CONTINUE
  200         CONTINUE
          ELSE
*
*           Form  C := alpha*A*B'          + beta*C
*
              DO 250 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 210 I = 1,M
                          C(I,J) = ZERO
  210                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 220 I = 1,M
                          C(I,J) = BETA*C(I,J)
  220                 CONTINUE
                  END IF
                  DO 240 L = 1,K
                      IF (B(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*B(J,L)
                          DO 230 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  230                     CONTINUE
                      END IF
  240             CONTINUE
  250         CONTINUE
          END IF
      ELSE IF (CONJA) THEN
          IF (CONJB) THEN
*
*           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
*
              DO 280 J = 1,N
                  DO 270 I = 1,M
                      TEMP = ZERO
                      DO 260 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*DCONJG(B(J,L))
  260                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  270             CONTINUE
  280         CONTINUE
          ELSE
*
*           Form  C := alpha*conjg( A' )*B' + beta*C
*
              DO 310 J = 1,N
                  DO 300 I = 1,M
                      TEMP = ZERO
                      DO 290 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*B(J,L)
  290                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  300             CONTINUE
  310         CONTINUE
          END IF
      ELSE
          IF (CONJB) THEN
*
*           Form  C := alpha*A'*conjg( B' ) + beta*C
*
              DO 340 J = 1,N
                  DO 330 I = 1,M
                      TEMP = ZERO
                      DO 320 L = 1,K
                          TEMP = TEMP + A(L,I)*DCONJG(B(J,L))
  320                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  330             CONTINUE
  340         CONTINUE
          ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
              DO 370 J = 1,N
                  DO 360 I = 1,M
                      TEMP = ZERO
                      DO 350 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  350                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  360             CONTINUE
  370         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of ZGEMM .
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX A(LDA,*),X(*),Y(*)
*     ..
*
*  Purpose
*  =======
*
*  ZGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
*
*     y := alpha*conjg( A' )*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE COMPLEX ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      DOUBLE COMPLEX ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*     .. Local Scalars ..
      DOUBLE COMPLEX TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
      LOGICAL NOCONJ
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +    .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGEMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
*
      NOCONJ = LSAME(TRANS,'T')
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      DO 50 I = 1,M
                          Y(I) = Y(I) + TEMP*A(I,J)
   50                 CONTINUE
                  END IF
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IY = KY
                      DO 70 I = 1,M
                          Y(IY) = Y(IY) + TEMP*A(I,J)
                          IY = IY + INCY
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
*
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 110 J = 1,N
                  TEMP = ZERO
                  IF (NOCONJ) THEN
                      DO 90 I = 1,M
                          TEMP = TEMP + A(I,J)*X(I)
   90                 CONTINUE
                  ELSE
                      DO 100 I = 1,M
                          TEMP = TEMP + DCONJG(A(I,J))*X(I)
  100                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  110         CONTINUE
          ELSE
              DO 140 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  IF (NOCONJ) THEN
                      DO 120 I = 1,M
                          TEMP = TEMP + A(I,J)*X(IX)
                          IX = IX + INCX
  120                 CONTINUE
                  ELSE
                      DO 130 I = 1,M
                          TEMP = TEMP + DCONJG(A(I,J))*X(IX)
                          IX = IX + INCX
  130                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  140         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of ZGEMV .
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX ALPHA
      INTEGER INCX,INCY,LDA,M,N
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX A(LDA,*),X(*),Y(*)
*     ..
*
*  Purpose
*  =======
*
*  ZGERU  performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Arguments
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE COMPLEX ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*     .. Local Scalars ..
      DOUBLE COMPLEX TEMP
      INTEGER I,INFO,IX,J,JY,KX
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGERU ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
          DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
      END IF
*
      RETURN
*
*     End of ZGERU .
*
      END
C
C
      DOUBLE PRECISION FUNCTION DCABS1(Z)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX Z
*     ..
*     ..
*  Purpose
*  =======
*
*  DCABS1 computes absolute value of a double complex number 
*
*     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DIMAG
*
      DCABS1 = ABS(DBLE(Z)) + ABS(DIMAG(Z))
      RETURN
      END
C
C
C======================================================================
C
      SUBROUTINE ZSCAL(N,ZA,ZX,INCX)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX ZA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX ZX(*)
*     ..
*
*  Purpose
*  =======
*
*     scales a vector by a constant.
*     jack dongarra, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,IX
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      IX = 1
      DO 10 I = 1,N
          ZX(IX) = ZA*ZX(IX)
          IX = IX + INCX
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
   20 DO 30 I = 1,N
          ZX(I) = ZA*ZX(I)
   30 CONTINUE
      RETURN
      END
C
C
C======================================================================
C
      SUBROUTINE ZSWAP(N,ZX,INCX,ZY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX ZX(*),ZY(*)
*     ..
*
*  Purpose
*  =======
*
*     interchanges two vectors.
*     jack dongarra, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      DOUBLE COMPLEX ZTEMP
      INTEGER I,IX,IY
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*       code for unequal increments or equal increments not equal
*         to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          ZTEMP = ZX(IX)
          ZX(IX) = ZY(IY)
          ZY(IY) = ZTEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*       code for both increments equal to 1
   20 DO 30 I = 1,N
          ZTEMP = ZX(I)
          ZX(I) = ZY(I)
          ZY(I) = ZTEMP
   30 CONTINUE
      RETURN
      END
C
C
C======================================================================
C
      SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX A(LDA,*),B(LDB,*)
*     ..
*
*  Purpose
*  =======
*
*  ZTRMM  performs one of the matrix-matrix operations
*
*     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
*
*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
*
*  Arguments
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE specifies whether  op( A ) multiplies B from
*           the left or right as follows:
*
*              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*
*              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain the matrix  B,  and  on exit  is overwritten  by the
*           transformed matrix.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
*     ..
*     .. Local Scalars ..
      DOUBLE COMPLEX TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      DOUBLE COMPLEX ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      DOUBLE COMPLEX ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*
*     Test the input parameters.
*
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOCONJ = LSAME(TRANSA,'T')
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRMM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (N.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*A*B.
*
              IF (UPPER) THEN
                  DO 50 J = 1,N
                      DO 40 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              DO 30 I = 1,K - 1
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   30                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP*A(K,K)
                              B(K,J) = TEMP
                          END IF
   40                 CONTINUE
   50             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              B(K,J) = TEMP
                              IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                              DO 60 I = K + 1,M
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   60                         CONTINUE
                          END IF
   70                 CONTINUE
   80             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*A'*B   or   B := alpha*conjg( A' )*B.
*
              IF (UPPER) THEN
                  DO 120 J = 1,N
                      DO 110 I = M,1,-1
                          TEMP = B(I,J)
                          IF (NOCONJ) THEN
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              DO 90 K = 1,I - 1
                                  TEMP = TEMP + A(K,I)*B(K,J)
   90                         CONTINUE
                          ELSE
                              IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                              DO 100 K = 1,I - 1
                                  TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  100                         CONTINUE
                          END IF
                          B(I,J) = ALPHA*TEMP
  110                 CONTINUE
  120             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = 1,M
                          TEMP = B(I,J)
                          IF (NOCONJ) THEN
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              DO 130 K = I + 1,M
                                  TEMP = TEMP + A(K,I)*B(K,J)
  130                         CONTINUE
                          ELSE
                              IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                              DO 140 K = I + 1,M
                                  TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  140                         CONTINUE
                          END IF
                          B(I,J) = ALPHA*TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*B*A.
*
              IF (UPPER) THEN
                  DO 200 J = N,1,-1
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 170 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  170                 CONTINUE
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
  200             CONTINUE
              ELSE
                  DO 240 J = 1,N
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 210 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  210                 CONTINUE
                      DO 230 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 220 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  220                         CONTINUE
                          END IF
  230                 CONTINUE
  240             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
*
              IF (UPPER) THEN
                  DO 280 K = 1,N
                      DO 260 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = ALPHA*A(J,K)
                              ELSE
                                  TEMP = ALPHA*DCONJG(A(J,K))
                              END IF
                              DO 250 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  250                         CONTINUE
                          END IF
  260                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = TEMP*A(K,K)
                          ELSE
                              TEMP = TEMP*DCONJG(A(K,K))
                          END IF
                      END IF
                      IF (TEMP.NE.ONE) THEN
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
  280             CONTINUE
              ELSE
                  DO 320 K = N,1,-1
                      DO 300 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = ALPHA*A(J,K)
                              ELSE
                                  TEMP = ALPHA*DCONJG(A(J,K))
                              END IF
                              DO 290 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  290                         CONTINUE
                          END IF
  300                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = TEMP*A(K,K)
                          ELSE
                              TEMP = TEMP*DCONJG(A(K,K))
                          END IF
                      END IF
                      IF (TEMP.NE.ONE) THEN
                          DO 310 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  310                     CONTINUE
                      END IF
  320             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRMM .
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX A(LDA,*),X(*)
*     ..
*
*  Purpose
*  =======
*
*  ZTRMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x,   or   x := conjg( A' )*x,
*
*  where x is an n element vector and  A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Arguments
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   x := A*x.
*
*              TRANS = 'T' or 't'   x := A'*x.
*
*              TRANS = 'C' or 'c'   x := conjg( A' )*x.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE COMPLEX ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*     .. Local Scalars ..
      DOUBLE COMPLEX TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOCONJ,NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (N.EQ.0) RETURN
*
      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x := A*x.
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 10 I = 1,J - 1
                              X(I) = X(I) + TEMP*A(I,J)
   10                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX
                  DO 40 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 30 I = 1,J - 1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX + INCX
   30                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX + INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 50 I = N,J + 1,-1
                              X(I) = X(I) + TEMP*A(I,J)
   50                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   60             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 70 I = N,J + 1,-1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX - INCX
   70                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX - INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
*
*        Form  x := A'*x  or  x := conjg( A' )*x.
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 110 J = N,1,-1
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 90 I = J - 1,1,-1
                              TEMP = TEMP + A(I,J)*X(I)
   90                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 100 I = J - 1,1,-1
                              TEMP = TEMP + DCONJG(A(I,J))*X(I)
  100                     CONTINUE
                      END IF
                      X(J) = TEMP
  110             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 140 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 120 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + A(I,J)*X(IX)
  120                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 130 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + DCONJG(A(I,J))*X(IX)
  130                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX - INCX
  140             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 170 J = 1,N
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 150 I = J + 1,N
                              TEMP = TEMP + A(I,J)*X(I)
  150                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 160 I = J + 1,N
                              TEMP = TEMP + DCONJG(A(I,J))*X(I)
  160                     CONTINUE
                      END IF
                      X(J) = TEMP
  170             CONTINUE
              ELSE
                  JX = KX
                  DO 200 J = 1,N
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 180 I = J + 1,N
                              IX = IX + INCX
                              TEMP = TEMP + A(I,J)*X(IX)
  180                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 190 I = J + 1,N
                              IX = IX + INCX
                              TEMP = TEMP + DCONJG(A(I,J))*X(IX)
  190                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX + INCX
  200             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRMV .
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX A(LDA,*),B(LDB,*)
*     ..
*
*  Purpose
*  =======
*
*  ZTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
*
*  The matrix X is overwritten on B.
*
*  Arguments
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
*     ..
*     .. Local Scalars ..
      DOUBLE COMPLEX TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      DOUBLE COMPLEX ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      DOUBLE COMPLEX ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*
*     Test the input parameters.
*
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOCONJ = LSAME(TRANSA,'T')
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRSM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (N.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*inv( A )*B.
*
              IF (UPPER) THEN
                  DO 60 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      END IF
                      DO 50 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      END IF
                      DO 90 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*inv( A' )*B
*           or    B := alpha*inv( conjg( A' ) )*B.
*
              IF (UPPER) THEN
                  DO 140 J = 1,N
                      DO 130 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          IF (NOCONJ) THEN
                              DO 110 K = 1,I - 1
                                  TEMP = TEMP - A(K,I)*B(K,J)
  110                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                          ELSE
                              DO 120 K = 1,I - 1
                                  TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
  120                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                          END IF
                          B(I,J) = TEMP
  130                 CONTINUE
  140             CONTINUE
              ELSE
                  DO 180 J = 1,N
                      DO 170 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          IF (NOCONJ) THEN
                              DO 150 K = I + 1,M
                                  TEMP = TEMP - A(K,I)*B(K,J)
  150                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                          ELSE
                              DO 160 K = I + 1,M
                                  TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
  160                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                          END IF
                          B(I,J) = TEMP
  170                 CONTINUE
  180             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*B*inv( A ).
*
              IF (UPPER) THEN
                  DO 230 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 190 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  190                     CONTINUE
                      END IF
                      DO 210 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              DO 200 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  200                         CONTINUE
                          END IF
  210                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 220 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  220                     CONTINUE
                      END IF
  230             CONTINUE
              ELSE
                  DO 280 J = N,1,-1
                      IF (ALPHA.NE.ONE) THEN
                          DO 240 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  240                     CONTINUE
                      END IF
                      DO 260 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              DO 250 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  250                         CONTINUE
                          END IF
  260                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 270 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  270                     CONTINUE
                      END IF
  280             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*inv( A' )
*           or    B := alpha*B*inv( conjg( A' ) ).
*
              IF (UPPER) THEN
                  DO 330 K = N,1,-1
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = ONE/A(K,K)
                          ELSE
                              TEMP = ONE/DCONJG(A(K,K))
                          END IF
                          DO 290 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  290                     CONTINUE
                      END IF
                      DO 310 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = A(J,K)
                              ELSE
                                  TEMP = DCONJG(A(J,K))
                              END IF
                              DO 300 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  300                         CONTINUE
                          END IF
  310                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 320 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  320                     CONTINUE
                      END IF
  330             CONTINUE
              ELSE
                  DO 380 K = 1,N
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = ONE/A(K,K)
                          ELSE
                              TEMP = ONE/DCONJG(A(K,K))
                          END IF
                          DO 340 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  340                     CONTINUE
                      END IF
                      DO 360 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = A(J,K)
                              ELSE
                                  TEMP = DCONJG(A(J,K))
                              END IF
                              DO 350 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  350                         CONTINUE
                          END IF
  360                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 370 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  370                     CONTINUE
                      END IF
  380             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRSM .
*
      END
C
C
C======================================================================
C
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      DLAMCH = RMACH
      FIRST  = .FALSE.
      RETURN
*
*     End of DLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      FIRST = .FALSE.
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
         FIRST = .FALSE.
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) INTEGER 
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5
*
      END
C
C
C======================================================================
C
      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*     
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
*
*  Purpose
*  =======
*
*       This program sets problem and machine dependent parameters
*       useful for xHSEQR and its subroutines. It is called whenever 
*       ILAENV is called with 12 <= ISPEC <= 16
*
*  Arguments
*  =========
*
*       ISPEC  (input) integer scalar
*              ISPEC specifies which tunable parameter IPARMQ should
*              return.
*
*              ISPEC=12: (INMIN)  Matrices of order nmin or less
*                        are sent directly to xLAHQR, the implicit
*                        double shift QR algorithm.  NMIN must be
*                        at least 11.
*
*              ISPEC=13: (INWIN)  Size of the deflation window.
*                        This is best set greater than or equal to
*                        the number of simultaneous shifts NS.
*                        Larger matrices benefit from larger deflation
*                        windows.
*
*              ISPEC=14: (INIBL) Determines when to stop nibbling and
*                        invest in an (expensive) multi-shift QR sweep.
*                        If the aggressive early deflation subroutine
*                        finds LD converged eigenvalues from an order
*                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
*                        then the next QR sweep is skipped and early
*                        deflation is applied immediately to the
*                        remaining active diagonal block.  Setting
*                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
*                        multi-shift QR sweep whenever early deflation
*                        finds a converged eigenvalue.  Setting
*                        IPARMQ(ISPEC=14) greater than or equal to 100
*                        prevents TTQRE from skipping a multi-shift
*                        QR sweep.
*
*              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
*                        a multi-shift QR iteration.
*
*              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
*                        following meanings.
*                        0:  During the multi-shift QR sweep,
*                            xLAQR5 does not accumulate reflections and
*                            does not use matrix-matrix multiply to
*                            update the far-from-diagonal matrix
*                            entries.
*                        1:  During the multi-shift QR sweep,
*                            xLAQR5 and/or xLAQRaccumulates reflections and uses
*                            matrix-matrix multiply to update the
*                            far-from-diagonal matrix entries.
*                        2:  During the multi-shift QR sweep.
*                            xLAQR5 accumulates reflections and takes
*                            advantage of 2-by-2 block structure during
*                            matrix-matrix multiplies.
*                        (If xTRMM is slower than xGEMM, then
*                        IPARMQ(ISPEC=16)=1 may be more efficient than
*                        IPARMQ(ISPEC=16)=2 despite the greater level of
*                        arithmetic work implied by the latter choice.)
*
*       NAME    (input) character string
*               Name of the calling subroutine
*
*       OPTS    (input) character string
*               This is a concatenation of the string arguments to
*               TTQRE.
*
*       N       (input) integer scalar
*               N is the order of the Hessenberg matrix H.
*
*       ILO     (input) INTEGER
*       IHI     (input) INTEGER
*               It is assumed that H is already upper triangular
*               in rows and columns 1:ILO-1 and IHI+1:N.
*
*       LWORK   (input) integer scalar
*               The amount of workspace available.
*
*  Further Details
*  ===============
*
*       Little is known about how best to choose these parameters.
*       It is possible to use different values of the parameters
*       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
*
*       It is probably best to choose different parameters for
*       different matrices and different parameters at different
*       times during the iteration, but this has not been
*       implemented --- yet.
*
*
*       The best choices of most of the parameters depend
*       in an ill-understood way on the relative execution
*       rate of xLAQR3 and xLAQR5 and on the nature of each
*       particular eigenvalue problem.  Experiment may be the
*       only practical way to determine which choices are most
*       effective.
*
*       Following is a list of default values supplied by IPARMQ.
*       These defaults may be adjusted in order to attain better
*       performance in any particular computational environment.
*
*       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
*                        Default: 75. (Must be at least 11.)
*
*       IPARMQ(ISPEC=13) Recommended deflation window size.
*                        This depends on ILO, IHI and NS, the
*                        number of simultaneous shifts returned
*                        by IPARMQ(ISPEC=15).  The default for
*                        (IHI-ILO+1).LE.500 is NS.  The default
*                        for (IHI-ILO+1).GT.500 is 3*NS/2.
*
*       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
*
*       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
*                        a multi-shift QR iteration.
*
*                        If IHI-ILO+1 is ...
*
*                        greater than      ...but less    ... the
*                        or equal to ...      than        default is
*
*                                0               30       NS =   2+
*                               30               60       NS =   4+
*                               60              150       NS =  10
*                              150              590       NS =  **
*                              590             3000       NS =  64
*                             3000             6000       NS = 128
*                             6000             infinity   NS = 256
*
*                    (+)  By default matrices of this order are
*                         passed to the implicit double shift routine
*                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
*                         values of NS are used only in case of a rare
*                         xLAHQR failure.
*
*                    (**) The asterisks (**) indicate an ad-hoc
*                         function increasing from 10 to 64.
*
*       IPARMQ(ISPEC=16) Select structured matrix multiply.
*                        (See ISPEC=16 above for details.)
*                        Default: 3.
*
*     ================================================================
*     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,
     $                   ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14,
     $                   NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
*     ..
*     .. Local Scalars ..
      INTEGER            NH, NS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
*     ..
*     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR.
     $    ( ISPEC.EQ.IACC22 ) ) THEN
*
*        ==== Set the number simultaneous shifts ====
*
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 )
     $      NS = 4
         IF( NH.GE.60 )
     $      NS = 10
         IF( NH.GE.150 )
     $      NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
         IF( NH.GE.590 )
     $      NS = 64
         IF( NH.GE.3000 )
     $      NS = 128
         IF( NH.GE.6000 )
     $      NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
*
      IF( ISPEC.EQ.INMIN ) THEN
*
*
*        ===== Matrices of order smaller than NMIN get sent
*        .     to xLAHQR, the classic double shift algorithm.
*        .     This must be at least 11. ====
*
         IPARMQ = NMIN
*
      ELSE IF( ISPEC.EQ.INIBL ) THEN
*
*        ==== INIBL: skip a multi-shift qr iteration and
*        .    whenever aggressive early deflation finds
*        .    at least (NIBBLE*(window size)/100) deflations. ====
*
         IPARMQ = NIBBLE
*
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
*
*        ==== NSHFTS: The number of simultaneous shifts =====
*
         IPARMQ = NS
*
      ELSE IF( ISPEC.EQ.INWIN ) THEN
*
*        ==== NW: deflation window size.  ====
*
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
*
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
*
*        ==== IACC22: Whether to accumulate reflections
*        .     before updating the far-from-diagonal elements
*        .     and whether to use 2-by-2 block structure while
*        .     doing it.  A small amount of work could be saved
*        .     by making this choice dependent also upon the
*        .     NH=IHI-ILO+1.
*
         IPARMQ = 0
         IF( NS.GE.KACMIN )
     $      IPARMQ = 1
         IF( NS.GE.K22MIN )
     $      IPARMQ = 2
*
      ELSE
*        ===== invalid value of ispec =====
         IPARMQ = -1
*
      END IF
*
*     ==== End of IPARMQ ====
*
      END
C
C
C======================================================================
C
      INTEGER FUNCTION IZAMAX(N,ZX,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX ZX(*)
*     ..
*
*  Purpose
*  =======
*
*     finds the index of element having max. absolute value.
*     jack dongarra, 1/15/85.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      DOUBLE PRECISION SMAX
      INTEGER I,IX
*     ..
*     .. External Functions ..
      DOUBLE PRECISION DCABS1
      EXTERNAL DCABS1
*     ..
      IZAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IZAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      IX = 1
      SMAX = DCABS1(ZX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          IF (DCABS1(ZX(IX)).LE.SMAX) GO TO 5
          IZAMAX = I
          SMAX = DCABS1(ZX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
   20 SMAX = DCABS1(ZX(1))
      DO 30 I = 2,N
          IF (DCABS1(ZX(I)).LE.SMAX) GO TO 30
          IZAMAX = I
          SMAX = DCABS1(ZX(I))
   30 CONTINUE
      RETURN
      END
C


