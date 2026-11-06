C
C=======================================================================
C
      SUBROUTINE SPEC_RAD_POWER(JE,E_KIN)
C
C   This subroutine computes the spectral radius of the G_o T kernel matrix
C       using the power method and the Rayleigh quotient. In addition,
C       a scalar convergence acceleration method can be used to refine the results.
C
C   Power method approximation of the spectral radius :
C
C      SR(A) = lim           ||A^{n} x|| / ||A^{n-1} x||
C              n --> +infty
C                                for any starting vector x
C
C   Rayleigh quotient approximation of the spectral radius :
C
C      SR(A) = lim           ||(A^{n}x)^{+} A  A^{n} x| / |(A^{n}x)^{+}  A^{n} x||
C              n --> +infty
C                                for any starting vector x
C
C   Note: a shift (= SHIFT) can be introduced in the power method.
C           In this case, only the Rayleigh quotient is valid
C
C        !!!!!  This version does not store the G_o T matrix !!!!!
C
C  Author : D. Sebilleau
C
C                                 Last modified :  1 Mar 2013
C
C      INCLUDE 'spec.inc'
      USE DIM_MOD
      USE CONVTYP_MOD
      USE COOR_MOD
      USE OUTFILES_MOD
      USE OUTUNITS_MOD
      USE TRANS_MOD
C
C!      PARAMETER(NLTWO=2*NL_M) !Moved to DIM_MOD
C
      INTEGER NN(0:N_ORD_M)
C
      REAL R_SN(N_ORD_M),R_SN1(N_ORD_M)
C
      COMPLEX*16 TLK,SUM_L1,SUM_L2,EXPKJ,EXPJK,ONEC,MULT
      COMPLEX*16 YLM(0:NLTWO,-NLTWO:NLTWO)
      COMPLEX*16 X(LINMAX*NATCLU_M)
      COMPLEX*16 AX(LINMAX*NATCLU_M),SUM_AX,XAX,XX
      COMPLEX*16 HL1(0:NLTWO)
      COMPLEX*16 ZEROC,IC,JC,JMULT,IFOURPI,CHIFT
      COMPLEX*16 SN(-1:N_ORD_M,-1:N_ORD_M)
      COMPLEX*16 COEF_GAUNT(LINMAX,LINMAX,0:N_GAUNT)
C
      REAL*8 PI,ATTKJ,GNT(0:N_GAUNT),XKJ,YKJ,ZKJ,RKJ,ZDKJ,KRKJ
      REAL*8 ATTJK,XJK,YJK,ZJK,RJK,ZDJK,KRJK,SMALL
      REAL*8 NORMX2,NORMAX2,RLIN
C
      CHARACTER*10 METH
      CHARACTER*24 OUTFILE,PATH
C
      DATA PI,SMALL /3.1415926535898D0,1.0D-6/
C
      ONEC=(1.D0,0.D0)
      IC=(0.0D0,1.0D0)
      ZEROC=(0.D0,0.D0)
      IFOURPI=4.D0*PI*IC
      IBESS=3
      N_CALL=0
      N_ACC=0
      K_ACC=0
      I_CONV=0
C
      CHIFT=DCMPLX(SHIFT)
C
      IF(N_MAX.GT.N_ORD_M) THEN
        WRITE(IUO1,99) N_MAX
        STOP
      ENDIF
C
C  Computing the size NLIN of the vectors A x
C
      LMAXJ=0
      NLIN=0
      DO JTYP=1,N_PROT
        NBTYPJ=NATYP(JTYP)
        LMJ=LMAX(JTYP,JE)
        LMAXJ=MAX(LMAXJ,LMJ)
        DO JNUM=1,NBTYPJ
          DO LJ=0,LMJ
            DO MJ=-LJ,LJ
              NLIN=NLIN+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C  Storage of the modified Gaunt coefficients
C
      LIN1=0
      DO L1=0,LMAXJ
        DO M1=-L1,L1
          LIN1=LIN1+1
          LIN2=0
          DO L2=0,LMAXJ
            L_MIN=ABS(L2-L1)
            L_MAX=L1+L2
            DO M2=-L2,L2
              LIN2=LIN2+1
              CALL GAUNT2(L1,M1,L2,M2,GNT)
              DO L=L_MIN,L_MAX,2
                M=M1-M2
                COEF_GAUNT(LIN1,LIN2,L)=IFOURPI*GNT(L)*(IC**
     &          L)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C  Starting vector for power method:
C
C               (JC^0,JC^1,JC^2,JC^3, ... ,JC^N)^EXPO
C
      JC=DCOS(PI/6.D0)+IC*DSIN(PI/6.D0)
      JMULT=JC**EXPO
C
      WRITE(IUO1,6)
C
      IF(ABS(SHIFT).LE.SMALL) THEN
        WRITE(36,76)
      ELSE
        WRITE(36,79)
      ENDIF
      WRITE(36,75)
C
C  Initialization the vectors and their squared norm
C
      MULT=JMULT
C
      AX(1)=ONEC
C
      NORMAX2=1.D0
C
      DO JLIN=2,NLIN
C
        AX(JLIN)=AX(JLIN-1)*MULT
        NORMAX2=NORMAX2+CDABS(AX(JLIN))*CDABS(AX(JLIN))
C
      ENDDO
C
C  Computation of the vectors and their squared norm
C
C            X ---> A*X
C
 120  IF(N_CALL.EQ.0) THEN
        J_INI=1
      ELSE
        J_INI=N_CALL*N_ITER+1
      ENDIF
      J_FIN=MIN((N_CALL+1)*N_ITER,N_MAX)
C
      IF(J_INI.GT.J_FIN) GOTO 112
C
      DO JORD=J_INI,J_FIN
C
        NORMX2=NORMAX2
C
        XX=ZEROC
        DO JLIN=1,NLIN
C
          X(JLIN)=AX(JLIN)
          XX=XX+CONJG(X(JLIN))*X(JLIN)
C
        ENDDO
C
        NORMAX2=0.D0
        XAX=ZEROC
C
C
C  Construction of the the multiple scattering kernel matrix G_o T
C                        times the vector x
C
C  Linear indices JLIN and KLIN represent (J,LJ) and (K,LK)
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
          LINJ=0
          DO LJ=0,LMJ
            DO MJ=-LJ,LJ
              JLIN=JLIN+1
              LINJ=LINJ+1
C
              SUM_AX=ZEROC
C
              KLIN=0
              DO KTYP=1,N_PROT
                NBTYPK=NATYP(KTYP)
                LMK=LMAX(KTYP,JE)
                DO KNUM=1,NBTYPK
                  KATL=NCORR(KNUM,KTYP)
                  XK=SYM_AT(1,KATL)
                  YK=SYM_AT(2,KATL)
                  ZK=SYM_AT(3,KATL)
C
                  IF(KATL.NE.JATL) THEN
C
                    XKJ=DBLE(XK-XJ)
                    YKJ=DBLE(YK-YJ)
                    ZKJ=DBLE(ZK-ZJ)
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
C
                  ENDIF
C
                  LINK=0
                  DO LK=0,LMK
                    L_MIN=ABS(LK-LJ)
                    L_MAX=LK+LJ
                    TLK=ATTKJ*DCMPLX(TL(LK,1,KTYP,JE))
                    DO MK=-LK,LK
                      KLIN=KLIN+1
                      LINK=LINK+1
                      IF(KATL.NE.JATL) THEN
C
                        SUM_L1=ZEROC
                        DO L=L_MIN,L_MAX,2
                          M=MJ-MK
                          IF(ABS(M).LE.L) THEN
                            SUM_L1=SUM_L1+HL1(L)*
     &                      COEF_GAUNT(LINJ,LINK,L)*YLM(L,M)
                          ENDIF
                        ENDDO
                        SUM_AX=SUM_AX+TLK*SUM_L1*X(
     &                  KLIN)
C
                      ENDIF
                      IF(KLIN.EQ.JLIN) SUM_AX=SUM_AX-
     &                CHIFT*X(KLIN)
C
                    ENDDO
                  ENDDO
C
                ENDDO
              ENDDO
C
              AX(JLIN)=SUM_AX
              NORMAX2=NORMAX2+CDABS(SUM_AX)*CDABS(SUM_AX)
              XAX=XAX+CONJG(X(JLIN))*SUM_AX
C
            ENDDO
          ENDDO
C
         ENDDO
       ENDDO
C
C  R_SN : power method
C  R_SN1: Rayleigh quotient
C
        R_SN(JORD)=SQRT(REAL(NORMAX2/NORMX2))
        R_SN1(JORD)=CDABS((XAX/XX)+SHIFT)
C
        IF(ABS(SHIFT).LE.SMALL) THEN
          WRITE(36,77) JORD,R_SN(JORD),R_SN1(JORD)
        ELSE
          WRITE(36,77) JORD,R_SN1(JORD)
        ENDIF
C
C  Check of convergence of power method series and
C          Rayleigh quotient method
C
C  Note: only Rayleigh quotient valid if shift
C
        IF((I_PWM.EQ.2).AND.(ABS(SHIFT).LE.SMALL)) THEN
          IF(JORD.GT.N_TABLE) THEN
            CALL CONV_SERIES(REAL(ACC),N_TABLE,JORD,R_SN,*555)
          ENDIF
        ENDIF
C
        IF((ABS(SHIFT).GT.SMALL).AND.(I_PWM.GT.0)) I_PWM=-I_PWM
        IF(I_PWM.EQ.-2) THEN
          IF(JORD.GT.N_TABLE) THEN
            CALL CONV_SERIES(REAL(ACC),N_TABLE,JORD,R_SN1,*555)
          ENDIF
        ENDIF
C
      ENDDO
C
  555 IF(J_INI.EQ.1) THEN
        WRITE(IUO1,10)
        WRITE(IUO1,11) SHIFT
        WRITE(IUO1,46) NLIN
        WRITE(IUO1,10)
        WRITE(IUO1,61)
      ENDIF
C
      WRITE(IUO1,10)
      IF(ABS(SHIFT).LE.SMALL) THEN
        WRITE(IUO1,44) R_SN(JORD-1)
      ENDIF
      WRITE(IUO1,53) R_SN1(JORD-1)
      WRITE(IUO1,47) JORD-1
      WRITE(IUO1,10)
      IF(J_FIN.EQ.N_MAX) THEN
        IF((I_PWM.EQ.1).OR.(I_PWM.EQ.3)) THEN
          WRITE(IUO1,213)
        ELSE
          WRITE(IUO1,60)
        ENDIF
      ENDIF
C
      IF(((I_ACC.EQ.1).OR.(I_ACC.EQ.2)).AND.(ABS(SHIFT).LE.SMALL)) THEN
C
C  Convergence acceleration on the power method results whenever necessary
C
        CALL ACC_CONV(N_CALL,J_INI,J_FIN,N_ACC,K_ACC,K_MAX,I_CONV,
     &  R_SN,SN,METH)
        IF(I_ACC.EQ.1) GOTO 111
        IF((I_CONV.EQ.1).AND.(I_ACC.EQ.2)) GOTO 111
        IF((I_CONV.EQ.0).AND.(I_ACC.EQ.2)) THEN
          IF(J_FIN.EQ.N_MAX) THEN
            GOTO 111
          ELSE
            GOTO 120
          ENDIF
        ENDIF
C
  111   WRITE(IUO1,10)
        WRITE(IUO1,61)
        WRITE(IUO1,10)
        WRITE(IUO1,52) METH
        WRITE(IUO1,48) CDABS(SN(N_ACC,K_ACC))
        WRITE(IUO1,49) N_ACC
        WRITE(IUO1,51) K_ACC
C
        IF(I_ACC.EQ.2) THEN
          OPEN(UNIT=35, FILE='div/n_k_table_pm.lis', STATUS=
     &    'unknown')
          K_MAX=MIN(K_MAX,J_FIN)
          DO L=0,K_MAX
            NN(L)=L
          ENDDO
          WRITE(35,*) '      '
          WRITE(35,78) E_KIN,METH
          WRITE(35,*) '      '
          WRITE(35,198) (NN(K),K=0,K_MAX)
          WRITE(35,199)
          WRITE(35,*) '      '
          DO N=0,K_MAX
            WRITE(35,200) N,(CDABS(SN(N,K)), K=0,K_MAX-N)
            WRITE(35,*) '      '
          ENDDO
        ENDIF
C
        WRITE(IUO1,10)
        IF(I_ACC.GE.1) WRITE(IUO1,60)
        IF(I_ACC.EQ.2) THEN
          IF(I_CONV.EQ.0) WRITE(IUO1,212)
          WRITE(IUO1,210)
        ENDIF
C
      ENDIF
C
      IF((ABS(SHIFT).GT.SMALL).AND.(I_ACC.GT.0)) I_ACC=-I_ACC
      IF((I_ACC.EQ.-1).OR.(I_ACC.EQ.-2))THEN
C
C  Convergence acceleration on the Rayleigh quotient whenever necessary
C
        CALL ACC_CONV(N_CALL,J_INI,J_FIN,N_ACC,K_ACC,K_MAX,I_CONV,
     &  R_SN1,SN,METH)
        IF(I_ACC.EQ.-1) GOTO 999
        IF((I_CONV.EQ.1).AND.(I_ACC.EQ.-2)) GOTO 999
        IF((I_CONV.EQ.0).AND.(I_ACC.EQ.-2)) THEN
          IF(J_FIN.EQ.N_MAX) THEN
            GOTO 999
          ELSE
            GOTO 120
          ENDIF
        ENDIF
C
  999   WRITE(IUO1,10)
        WRITE(IUO1,61)
        WRITE(IUO1,10)
        WRITE(IUO1,52) METH
        WRITE(IUO1,48) CDABS(SN(N_ACC,K_ACC))
        WRITE(IUO1,49) N_ACC
        WRITE(IUO1,51) K_ACC
C
        IF(I_ACC.EQ.-2) THEN
          OPEN(UNIT=35, FILE='div/n_k_table_rq.lis', STATUS=
     &    'unknown')
          K_MAX=MIN(K_MAX,J_FIN)
          DO L=0,K_MAX
            NN(L)=L
          ENDDO
          WRITE(35,*) '      '
          WRITE(35,78) E_KIN,METH
          WRITE(35,*) '      '
          WRITE(35,198) (NN(K),K=0,K_MAX)
          WRITE(35,199)
          WRITE(35,*) '      '
          DO N=0,K_MAX
            WRITE(35,200) N,(CDABS(SN(N,K)), K=0,K_MAX-N)
            WRITE(35,*) '      '
          ENDDO
        ENDIF
C
        WRITE(IUO1,10)
        IF(ABS(I_ACC).GE.1) WRITE(IUO1,60)
        IF(I_ACC.EQ.-2) THEN
          IF(I_CONV.EQ.0) WRITE(IUO1,215)
          WRITE(IUO1,214)
        ENDIF
C
      ENDIF
C
C  Writing the results:
C
      IF(ABS(I_ACC).LT.2) THEN
C
C  No convergence acceleration: power method and Rayleigh quotient
C
        IF(ABS(SHIFT).LE.SMALL) THEN
          WRITE(IUO2,*) E_KIN,R_SN(JORD-1),R_SN1(JORD-1)
        ELSE
          WRITE(IUO2,*) E_KIN,R_SN1(JORD-1)
        ENDIF
C
      ELSE
C
        IF(I_CONV.EQ.1) THEN
C
C  Convergence acceleration has succeeded:
C
          WRITE(IUO2,*) E_KIN,REAL(CDABS(SN(N_ACC,K_ACC)))
C
        ELSEIF(I_CONV.EQ.0) THEN
C
C  Convergence acceleration has failed: power method and
C         Rayleigh quotient used instead
C
          IF(ABS(SHIFT).LE.SMALL) THEN
            WRITE(IUO2,*) E_KIN,R_SN(JORD-1),R_SN1(JORD-1)
          ELSE
            WRITE(IUO2,*) E_KIN,R_SN1(JORD-1)
          ENDIF
C
        ENDIF
      ENDIF
C
 112  WRITE(IUO1,211)
C
      RETURN
C
   5  FORMAT(/,11X,'-----------------  EIGENVALUE ANALYSIS ','---------
     &--------')
   6  FORMAT(/,11X,'-----------  SHIFTED POWER METHOD ANALYSIS ',' ----
     &--------')
  10  FORMAT(11X,'-',54X,'-')
  11  FORMAT(11X,'-',4X,'SHIFT            ',20X,':   ',F6.3,3X,'-')
  15  FORMAT(11X,'-',14X,'MAXIMUM MODULUS : ',F9.6,13X,'-')
  20  FORMAT(11X,'-',14X,'MINIMUM MODULUS : ',F9.6,13X,'-')
  25  FORMAT(11X,'-',6X,'1 EIGENVALUE IS > 1 ON A TOTAL OF ',I8,6X,'-')
  30  FORMAT(11X,'-',4X,I5,' EIGENVALUES ARE > 1 ON A TOTAL OF ',I8,2X,
     &'-')
  35  FORMAT(11X,'-',11X,'THE ',I3,' LARGER EIGENVALUES ARE :',11X,'-')
  40  FORMAT(11X,'-',6X,F7.4,2X,F7.4,2X,F7.4,2X,F7.4,2X,F7.4,5X,'-')
  44  FORMAT(11X,'-',4X,'SPECTRAL RADIUS BY THE POWER METHOD  : ',F8.5,
     &3X,'-')
  45  FORMAT(11X,'-',4X,'SPECTRAL RADIUS OF THE KERNEL MATRIX : ',F8.5,
     &3X,'-')
  46  FORMAT(11X,'-',4X,'MATRIX SIZE                          : ',I8,
     &3X,'-')
  47  FORMAT(11X,'-',4X,'POWER METHOD TRUNCATION ORDER        : ',I8,
     &3X,'-')
  48  FORMAT(11X,'-',4X,'SPECTRAL RADIUS BY ACCELERATION      : ',F8.5,
     &3X,'-')
  49  FORMAT(11X,'-',4X,'ACCELERATION TRUNCATION ORDER N      : ',I8,
     &3X,'-')
  50  FORMAT(11X,'-',4X,'---> THE MULTIPLE SCATTERING SERIES ',
     &'CONVERGES ',4X,'-')
  51  FORMAT(11X,'-',4X,'ACCELERATION TRUNCATION ORDER K      : ',I8,
     &3X,'-')
  52  FORMAT(11X,'-',4X,'ACCELERATION METHOD                  : ',A10,
     &1X,'-')
  53  FORMAT(11X,'-',4X,'SPECTRAL RADIUS BY RAYLEIGH QUOTIENT : ',F8.5,
     &3X,'-')
  55  FORMAT(11X,'-',10X,'---> NO CONVERGENCE OF THE MULTIPLE',9X,'-',/
     &,11X,'-',18X,'SCATTERING SERIES',19X,'-')
  60  FORMAT(11X,'----------------------------------------','----------
     &------',/)
  61  FORMAT(11X,'----------------------------------------','----------
     &------')
  65  FORMAT(11X,'-',4X,'    LABEL OF LARGEST EIGENVALUE   :  ',I5,8X,
     &'-')
  70  FORMAT(11X,'-',4X,'    LARGEST EIGENVALUE   : ','(',F6.3,',',F6.
     &3,')',8X,'-')
  75  FORMAT('              ')
  76  FORMAT(2X,'ITERATION:      SPECTRAL RADIUS:','      RAYLEIGH 
     &QUOTIENT:')
  77  FORMAT(2X,I6,15X,F6.3,15X,F6.3)
  78  FORMAT(10X,'    E_KIN = ',F6.2,' eV  -- (N,K) TABLE OF ',A10,' 
     &METHOD')
  79  FORMAT(2X,'ITERATION:      RAYLEIGH QUOTIENT:')
  80  FORMAT('    KINETIC ENERGY : ',F7.2,' eV')
  85  FORMAT('    LARGEST MODULUS OF EIGENVALUE : ',F6.3)
  90  FORMAT('    LABEL OF LARGEST EIGENVALUE   : ',I5)
  95  FORMAT('    LARGEST EIGENVALUE            :  (',F6.3,',',F6.3,')
     &')
  99  FORMAT(///,12X,' <<<<<<<<<<  DIMENSION OF N_ORD_M IN INCLUDE',' 
     &FILE  >>>>>>>>>>',/,12X,' <<<<<<<<<<        SHOULD BE AT LEAST ',
     &I5,  6X,'  >>>>>>>>>>',///)
 100  FORMAT(5X,F7.3,2X,F7.3,2X,F6.3)
 105  FORMAT(7X,'EIGENVALUES :',3X,'MODULUS :')
 110  FORMAT(2X,'-------------------------------')
 198  FORMAT('    K',50(I3,4X))
 199  FORMAT('  N')
 200  FORMAT(I3,50(2X,F5.3))
 210  FORMAT(//,'   --->  THE (N,K) TABLE HAS BEEN WRITTEN ','IN /div/
     &n_k_table_pm.lis',//)
 211  FORMAT(//,'   --->  THE EVOLUTION OF THE SPECTRAL RADIUS ','HAS 
     &BEEN WRITTEN IN /div/radius_conv.lis',//)
 212  FORMAT(//,10X,'   !!!!! THE POWER METHOD HAS NOT CONVERGED YET  
     &','!!!!!')
 213  FORMAT(11X,'------------  NO CONVERGENCE CHECK WAS MADE ','------
     &------')
 214  FORMAT(//,'   --->  THE (N,K) TABLE HAS BEEN WRITTEN ','IN /div/
     &n_k_table_rq.lis',//)
 215  FORMAT(//,8X,'   !!!!! THE RAYLEIGH QUOTIENT HAS NOT CONVERGED',
     &' YET  !!!!!')
C
      END
