C
C=======================================================================
C
      SUBROUTINE SPEC_RAD_POWER(JE,E_KIN)
C
C   This subroutine stores the G_o T kernel matrix and computes
C       its spectral radius using the power method. In addition,
C       a scalar convergence acceleration method can be used to
C       refine the results.
C
C  Author : D. Sebilleau
C
C                                 Last modified :  4 Feb 2013
C
C      INCLUDE 'spec.inc'
      USE DIM_MOD
      USE CONVTYP_MOD
      USE COOR_MOD
      USE OUTFILES_MOD
      USE OUTUNITS_MOD
      USE TRANS_MOD
C
C
C!      PARAMETER(NLTWO=2*NL_M) !Moved to DIM_MOD
C
      INTEGER NN(0:N_ORD_M)
C
      REAL R_SN(N_ORD_M),R_SN1(N_ORD_M),R_SN2(N_ORD_M)
C
      COMPLEX*16 TLK,SUM_L,EXPKJ,ONEC,MULT
      COMPLEX*16 YLM(0:NLTWO,-NLTWO:NLTWO)
      COMPLEX*16 X(LINMAX*NATCLU_M),Y(LINMAX*NATCLU_M)
      COMPLEX*16 AX(LINMAX*NATCLU_M),AY(LINMAX*NATCLU_M),SUM_AX,SUM_AY
      COMPLEX*16 HL1(0:NLTWO),A(LINMAX*NATCLU_M,LINMAX*NATCLU_M)
      COMPLEX*16 ZEROC,IC,JC,JMULT
      COMPLEX*16 SN(-1:N_ORD_M,-1:N_ORD_M)
C
      REAL*8 PI,ATTKJ,GNT(0:N_GAUNT),XKJ,YKJ,ZKJ,RKJ,ZDKJ,KRKJ
      REAL*8 NORMX2,NORMY2,NORMAX2,NORMAY2,RLIN
C
C
      CHARACTER*10 METH
      CHARACTER*24 OUTFILE,PATH
C
      DATA PI /3.1415926535898D0/
C
      ONEC=(1.D0,0.D0)
      IC=(0.0D0,1.0D0)
      ZEROC=(0.D0,0.D0)
      IBESS=3
      N_CALL=0
      N_ACC=0
      K_ACC=0
      I_CONV=0
C
      IF(N_MAX.GT.N_ORD_M) THEN
        WRITE(IUO1,99) N_MAX
        STOP
      ENDIF
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
C  Construction of the the multiple scattering kernel matrix G_o T
C
C       Elements are stored in A(KLIN,JLIN) using linear indices
C          JLIN and KLIN  representing (J,LJ) and (K,LK)
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
                      A(KLIN,JLIN)=ZEROC
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
                      IF(KATL.NE.JATL) THEN
                        A(KLIN,JLIN)=TLK*SUM_L
                      ENDIF
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
      NLIN=JLIN
C
C  Power method approximation of the spectral radius :
C
C      SR(A) = lim           ||A x||^n / ||x||^n
C              n --> +infty
C                                for any starting vector x
C
      RLIN=DFLOAT(NLIN)
C
C  Initialization the vectors and their squared norm
C
      MULT=JMULT
C
      AX(1)=ONEC
      AY(1)=ONEC
C
      NORMAX2=1.D0
      NORMAY2=1.D0
C
      DO ILIN=2,NLIN
C
        AX(ILIN)=AX(ILIN-1)*MULT
        AY(ILIN)=AY(ILIN-1)*MULT
        NORMAX2=NORMAX2+DREAL(AX(ILIN)*DCONJG(AX(ILIN)))
        NORMAY2=NORMAY2+DREAL(AY(ILIN)*DCONJG(AY(ILIN)))
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
        NORMY2=NORMAY2
C
        DO JLIN=1,NLIN
C
          X(JLIN)=AX(JLIN)
          Y(JLIN)=AY(JLIN)
C
        ENDDO
C
        NORMAX2=0.D0
        NORMAY2=0.D0
C
        DO ILIN=1,NLIN
C
          SUM_AX=ZEROC
          SUM_AY=ZEROC
C
          DO JLIN=1,NLIN
C
            SUM_AX=SUM_AX+A(ILIN,JLIN)*X(JLIN)
            SUM_AY=SUM_AY+A(JLIN,ILIN)*Y(JLIN)
C
          ENDDO
C
          AX(ILIN)=SUM_AX
          AY(ILIN)=SUM_AY
          NORMAX2=NORMAX2+DREAL(SUM_AX*DCONJG(SUM_AX))
          NORMAY2=NORMAY2+DREAL(SUM_AY*DCONJG(SUM_AY))
C
        ENDDO
C
        R_SN1(JORD)=SQRT(REAL(NORMAX2/NORMX2))
        R_SN2(JORD)=SQRT(REAL(NORMAY2/NORMY2))
        R_SN(JORD)=SQRT(R_SN1(JORD)*R_SN2(JORD))
C
      ENDDO
C
      IF(J_INI.EQ.1) THEN
        WRITE(IUO1,10)
        WRITE(IUO1,46) NLIN
        WRITE(IUO1,10)
        WRITE(IUO1,61)
      ENDIF
C
      WRITE(IUO1,10)
      WRITE(IUO1,44) R_SN(JORD-1)
      WRITE(IUO1,47) JORD-1
C
  112 IF(I_ACC.GE.1) THEN
C
C  Convergence acceleration on the results whenever necessary
C
        CALL ACC_CONV(N_CALL,J_INI,J_FIN,N_ACC,K_ACC,K_MAX,I_CONV,
     &  R_SN,SN,METH)
        IF((I_CONV.EQ.1).AND.(I_ACC.EQ.1)) GOTO 111
        IF((I_CONV.EQ.0).OR.(I_ACC.EQ.2)) GOTO 120
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
          OPEN(UNIT=35, FILE='div/n_k_table.lis', STATUS='unknown')
          DO L=0,K_MAX
            NN(L)=L
          ENDDO
          WRITE(35,*) '      '
          WRITE(35,*) '                   (N,K) TABLE OF ',METH,' 
     &    METHOD'
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
        WRITE(IUO1,60)
        IF(I_ACC.EQ.2) THEN
          WRITE(IUO1,210)
        ENDIF
C
      ENDIF
C
      WRITE(IUO2,*) E_KIN,R_SN(JORD-1)
C
      RETURN
C
   5  FORMAT(/,11X,'-----------------  EIGENVALUE ANALYSIS ','---------
     &--------')
   6  FORMAT(/,11X,'----------------  POWER METHOD ANALYSIS ','--------
     &--------')
  10  FORMAT(11X,'-',54X,'-')
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
 210  FORMAT(//,'          --->  THE (N,K) TABLE HAS BEEN WRITTEN ',
     &'IN /div/n_k_table.lis',//)
C
      END
