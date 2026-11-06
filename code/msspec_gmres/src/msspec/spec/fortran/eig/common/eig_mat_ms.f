C
C=======================================================================
C
      SUBROUTINE EIG_MAT_MS(JE,E_KIN)
C
C   This subroutine stores the G_o T kernel matrix and computes
C     its eigenvalues to check the larger ones. If one or more
C     of these eigenvalues are larger than or equal to 1, no
C     series expansion can converge
C
C                                       Last modified : 20 Jul 2009
C
C      INCLUDE 'spec.inc'
C
      USE DIM_MOD
      USE COOR_MOD
      USE OUTFILES_MOD
      USE OUTUNITS_MOD
      USE TRANS_MOD
CKMD
      USE RENORM_MOD
CKMD
C
C!      PARAMETER(NLTWO=2*NL_M) !Moved to DIM_MOD
C
      CHARACTER*24 OUTFILE,PATH
      CHARACTER*512 MATRIX_FILE
      INTEGER MATRIX_IO,MATRIX_STATUS,MATRIX_COL
C
      COMPLEX*16 HL1(0:NLTWO),SM(LINMAX*NATCLU_M,LINMAX*NATCLU_M)
CKMD      COMPLEX*16 SUM_L,IC,ZEROC,WORK(32*LINMAX*NATCLU_M)
      COMPLEX*16 SUM_L,IC,ZEROC
      COMPLEX*16 YLM(0:NLTWO,-NLTWO:NLTWO),TLK,EXPKJ
      COMPLEX*16 W(LINMAX*NATCLU_M)
CKMD      COMPLEX*16 VL(1,1),VR(1,1)
C
CKMD      DOUBLE PRECISION RWORK(2*LINMAX*NATCLU_M)
C
      REAL*8 PI,ATTKJ,GNT(0:N_GAUNT),XKJ,YKJ,ZKJ,RKJ,ZDKJ,KRKJ
C
      REAL W1(LINMAX*NATCLU_M),W2(LINMAX*NATCLU_M)
C
      DATA PI,SMALL /3.1415926535898D0,0.0001/
C
      IC=(0.D0,1.D0)
      ZEROC=(0.D0,0.D0)
      IBESS=3
      NPRINT=10
C
      WRITE(IUO1,5)
CKMD
      IOUT2=55
      IOUT3=60
CKMD
      IF(JE.EQ.1) THEN
C
C  Name of the second output file where all the eigenvalues
C    will be written
C
CKMD Moved initialisation of IOUT2 and IOUT3 to before the test on JE,
CKMD otherwise they were not initialised for JE > 1
CKMD        IOUT2=55
CKMD        IOUT3=60
        N_DOT=1
        DO J_CHAR=1,24
          IF(OUTFILE2(J_CHAR:J_CHAR).EQ.'.') GOTO 888
          N_DOT=N_DOT+1
        ENDDO
  888   CONTINUE
        OUTFILE=OUTFILE2(1:N_DOT)//'egv'
        PATH=OUTFILE2(1:N_DOT)//'pth'
        OPEN(UNIT=IOUT2, FILE=OUTFILE, STATUS='UNKNOWN')
        OPEN(UNIT=IOUT3, FILE=PATH, STATUS='UNKNOWN')
      ENDIF
C
C  Construction of the multiple scattering kernel matrix G_o T.
C       Elements are stored using a linear index LINJ
C                    representing (J,LJ)
CKMD
      SM = CMPLX(0.0D0, 0.0D0)
C
      JLIN=0
      DO JTYP=1,N_PROT
        NBTYPJ=NATYP(JTYP)
        LMJ=LMAX(JTYP,JE)
        DO JNUM=1,NBTYPJ
          JATL=NCORR(JNUM,JTYP)
          XJ=SYM_AT(1,JATL)
          YJ=SYM_AT(2,JATL)
          ZJ=SYM_AT(3,JATL)
C
          DO LJ=0,LMJ
            DO MJ=-LJ,LJ
              JLIN=JLIN+1
C
              KLIN=0
              DO KTYP=1,N_PROT
                NBTYPK=NATYP(KTYP)
                LMK=LMAX(KTYP,JE)
                DO KNUM=1,NBTYPK
                  KATL=NCORR(KNUM,KTYP)
                  IF(KATL.NE.JATL) THEN
                    XKJ=DBLE(SYM_AT(1,KATL)-XJ)
                    YKJ=DBLE(SYM_AT(2,KATL)-YJ)
                    ZKJ=DBLE(SYM_AT(3,KATL)-ZJ)
                    RKJ=DSQRT(XKJ*XKJ+YKJ*YKJ+ZKJ*ZKJ)
                    KRKJ=DBLE(VK(JE))*RKJ
                    ATTKJ=DEXP(-DIMAG(DCMPLX(VK(JE)))*
     &              RKJ)
                    EXPKJ=(XKJ+IC*YKJ)/RKJ
                    ZDKJ=ZKJ/RKJ
                    CALL SPH_HAR2(2*NL_M,ZDKJ,EXPKJ,YLM,
     &              LMJ+LMK)
                    CALL BESPHE2(LMJ+LMK+1,IBESS,KRKJ,
     &              HL1)
                  ENDIF
C
                  DO LK=0,LMK
                    L_MIN=ABS(LK-LJ)
                    L_MAX=LK+LJ
                    TLK=DCMPLX(TL(LK,1,KTYP,JE))
                    DO MK=-LK,LK
                      KLIN=KLIN+1
                      SM(KLIN,JLIN)=ZEROC
                      SUM_L=ZEROC
                      IF(KATL.NE.JATL) THEN
                        CALL GAUNT2(LK,MK,LJ,MJ,GNT)
C
                        DO L=L_MIN,L_MAX,2
                          M=MJ-MK
                          IF(ABS(M).LE.L) THEN
                            SUM_L=SUM_L+(IC**L)*
     &                      HL1(L)*YLM(L,M)*GNT(L)
                          ENDIF
                        ENDDO
                        SUM_L=SUM_L*ATTKJ*4.D0*PI*IC
                      ELSE
                        SUM_L=ZEROC
                      ENDIF
C
                      SM(KLIN,JLIN)=TLK*SUM_L
C
                    ENDDO
                  ENDDO
C
                ENDDO
              ENDDO
C
            ENDDO
          ENDDO
C
        ENDDO
      ENDDO
C
      N_DIM=LINMAX*NATCLU_M
C
CKMD Optional diagnostic export for iterative eigensolver prototypes.
CKMD Disabled unless MSSPEC_KERNEL_MATRIX_FILE names an output path.
C
      MATRIX_FILE=' '
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_KERNEL_MATRIX_FILE',
     &     MATRIX_FILE,STATUS=MATRIX_STATUS)
      IF (MATRIX_STATUS.EQ.0.AND.LEN_TRIM(MATRIX_FILE).GT.0) THEN
        OPEN(NEWUNIT=MATRIX_IO,FILE=TRIM(MATRIX_FILE),STATUS='REPLACE',
     &       FORM='UNFORMATTED',ACCESS='STREAM',ACTION='WRITE',
     &       IOSTAT=MATRIX_STATUS)
        IF (MATRIX_STATUS.EQ.0) THEN
          WRITE(MATRIX_IO) JLIN
          DO MATRIX_COL=1,JLIN
            WRITE(MATRIX_IO) SM(1:JLIN,MATRIX_COL)
          ENDDO
          CLOSE(MATRIX_IO)
        ENDIF
      ENDIF
C
      IF (I_REN.gt.0) THEN
C
CKMD Renormalize the matrix SM
C
          CALL RENORM_MATRIX(JLIN,SM,N_DIM)
C
CKMD SM now contains the renormalized matrix
C
      END IF
C
C   Eigenvalues of the kernel multiple scattering matrix SM
C
CKMDC      CALL ZGEEV('N','N',JLIN,SM,N_DIM,W,VL,1,VR,1,WORK,32*N_DIM,RWORK,
CKMD          LWORK = 32*N_DIM
CKMD          CALL ZGEEV('N','N',JLIN,SM,N_DIM,W,VL,1,VR,1,WORK,LWORK,
CKMD     &RWORK,INFO)
CKMD
      CALL DIAG_MAT(JLIN,SM,N_DIM,W,INFO)
CKMD
CKMD SM has been overwritten here
C
CKMD      IF(INFO.NE.0) THEN
CKMD          WRITE(IUO1,*) '                        '
CKMD          WRITE(IUO1,*) '    --->  WORK(1),INFO =',WORK(1),INFO
CKMD          WRITE(IUO1,*) '                        '
CKMD      ENDIF
C
CKMD Save eigenvalues to unformatted stream file eigenvalues.dat
C
      call save_eigenvalues(w, jlin, e_kin)
C
      N_EIG=0
C
      WRITE(IOUT2,75)
      WRITE(IOUT2,110)
      WRITE(IOUT2,80) E_KIN
      WRITE(IOUT2,110)
      WRITE(IOUT2,75)
      WRITE(IOUT2,105)
      WRITE(IOUT2,75)
      XMAX_L=0.
      N_XMAX=0
C
      DO LIN=1,JLIN
        EIG=REAL(CDABS(W(LIN)))
        WRITE(IOUT2,100) REAL(DBLE(W(LIN))),REAL(DIMAG(W(LIN))),EIG
        IF((EIG-XMAX_L).GT.0.0001) N_XMAX=LIN
        XMAX_L=MAX(XMAX_L,EIG)
        W1(LIN)=EIG
        IF(EIG.GT.1.000) THEN
          N_EIG=N_EIG+1
        ENDIF
      ENDDO
C
      WRITE(IOUT2,75)
      WRITE(IOUT2,85) XMAX_L
      WRITE(IOUT2,90) N_XMAX
      WRITE(IOUT2,95) W(N_XMAX)
      WRITE(IOUT2,75)
C
      CALL ORDRE(JLIN,W1,NFIN,W2)
C
C
CKMD
      WRITE(IUO1,10)
      WRITE(IUO1,12) JLIN
CKMD
      WRITE(IUO1,10)
      WRITE(IUO1,10)
      WRITE(IUO1,15) W2(1)
      WRITE(IUO1,20) W2(NFIN)
      WRITE(IUO1,10)
      WRITE(IUO1,10)
      IF(N_EIG.GE.1) THEN
        IF(N_EIG.EQ.1) THEN
          WRITE(IUO1,25) JLIN
        ELSE
          WRITE(IUO1,30) N_EIG,JLIN
        ENDIF
      ENDIF
C
      WRITE(IUO1,65) N_XMAX
      WRITE(IUO1,70) W(N_XMAX)
      WRITE(IUO1,10)
      WRITE(IOUT3,100) REAL(DBLE(W(N_XMAX))),REAL(DIMAG(W(N_XMAX)))
C
      NPR=NPRINT/5
      WRITE(IUO1,10)
      WRITE(IUO1,10)
      WRITE(IUO1,35) 5*NPR
      WRITE(IUO1,10)
C
      DO JP=0,NPR-1
        J=5*JP
        WRITE(IUO1,40) W2(J+1),W2(J+2),W2(J+3),W2(J+4),W2(J+5)
      ENDDO
      WRITE(IUO1,10)
      WRITE(IUO1,10)
CKMD
      IF (I_REN.NE.0) THEN
         WRITE(IUO1,46) REN_R, REN_I
         WRITE(IUO1,47) W2(1)
      ELSE
         WRITE(IUO1,45) W2(1)
      ENDIF
CKMD
      IF (I_REN.NE.0) THEN
         WRITE(IUO2,*) E_KIN,W2(1) !,REN_R,REN_I
      ELSE
         WRITE(IUO2,*) E_KIN,W2(1)
      ENDIF
      IF(N_EIG.EQ.0) THEN
        WRITE(IUO1,50)
      ELSE
        WRITE(IUO1,55)
      ENDIF
      WRITE(IUO1,10)
      WRITE(IUO1,10)
      WRITE(IUO1,60)
C
      RETURN
C
CKMD
   5  FORMAT(/,11X,'-----------------  EIGENVALUE ANALYSIS ',
     &'-----------------')
  10  FORMAT(11X,'-',54X,'-')
CKMD
  12  FORMAT(11X,'-',14X,'MATRIX DIMENSION : ',I8,13X,'-')
CKMD
  15  FORMAT(11X,'-',14X,'MAXIMUM MODULUS : ',F9.6,13X,'-')
  20  FORMAT(11X,'-',14X,'MINIMUM MODULUS : ',F9.6,13X,'-')
  25  FORMAT(11X,'-',6X,'1 EIGENVALUE IS > 1 ON A TOTAL OF ',I8,6X,'-')
  30  FORMAT(11X,'-',4X,I5,' EIGENVALUES ARE > 1 ON A TOTAL OF ',I8,2X,
     &'-')
  35  FORMAT(11X,'-',11X,'THE ',I3,' LARGER EIGENVALUES ARE :',11X,'-')
  40  FORMAT(11X,'-',6X,F7.4,2X,F7.4,2X,F7.4,2X,F7.4,2X,F7.4,5X,'-')
  45  FORMAT(11X,'-',5X,'SPECTRAL RADIUS OF THE KERNEL MATRIX :',F6.3,
     &5X,'-')
CKMD
  46  FORMAT(11X,'-',16X,'OMEGA = (',F6.3,',',F6.3,')',15X,'-')
  47  FORMAT(11X,'-',2X,'SPECTRAL RADIUS OF THE RENORMALIZED MATRIX: ',
     &F6.3,2X,'-')
CKMD
  50  FORMAT(11X,'-',5X,'---> THE MULTIPLE SCATTERING SERIES ',
     &'CONVERGES',4X,'-')
  55  FORMAT(11X,'-',10X,'---> NO CONVERGENCE OF THE MULTIPLE',9X,'-',/
     &,11X,'-',18X,'SCATTERING SERIES',19X,'-')
CKMD
  60  FORMAT(11X,'----------------------------------------',
     &'----------------',/)
  65  FORMAT(11X,'-',5X,'    LABEL OF LARGEST EIGENVALUE   : ',I5,8X,'-
     &')
  70  FORMAT(11X,'-',5X,'    LARGEST EIGENVALUE  : ','(',F6.3,',',F6.3,
     &')',8X,'-')
  75  FORMAT('              ')
  80  FORMAT('    KINETIC ENERGY : ',F7.2,' eV')
  85  FORMAT('    LARGEST MODULUS OF EIGENVALUE : ',F6.3)
  90  FORMAT('    LABEL OF LARGEST EIGENVALUE   : ',I5)
  95  FORMAT('    LARGEST EIGENVALUE            :  (',F6.3,',',F6.3,')
     &')
 100  FORMAT(5X,F7.3,2X,F7.3,2X,F6.3)
 105  FORMAT(7X,'EIGENVALUES :',3X,'MODULUS :')
 110  FORMAT(2X,'-------------------------------')
C
      END
