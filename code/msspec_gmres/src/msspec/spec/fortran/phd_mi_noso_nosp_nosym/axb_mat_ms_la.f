C
C=======================================================================
C
      SUBROUTINE INV_MAT_MS(JE,TAU,SRC_COMP,NCOMP)
C
C   This subroutine stores the multiple scattering matrix and computes 
C     the scattering path operator TAU^{j 0} exactly, without explicitely 
C     using the inverse matrix.
C
C                           (Photoelectron case)
C
C                                         Last modified : 28 Mar 2007
C
      USE DIM_MOD
C
      USE COOR_MOD
      USE INIT_L_MOD
      USE TRANS_MOD

C      PARAMETER(NLTWO=2*NL_M)
C
      COMPLEX*16 HL1(0:2*NL_M)
      COMPLEX*16, ALLOCATABLE :: SM(:,:)
      COMPLEX*16 IN(LINMAX*NATCLU_M,LINMAX)
      COMPLEX SRC_COMP(LINFMAX,LINFMAX)
      COMPLEX*16 SUM_L,ONEC,IC,ZEROC
      COMPLEX*16 YLM(0:2*NL_M,-2*NL_M:2*NL_M),TLJ,TLK,EXPKJ
C
      COMPLEX TAU(LINMAX,LINFMAX,NATCLU_M)
C
      REAL*8 PI,ATTKJ,GNT(0:N_GAUNT),XKJ,YKJ,ZKJ,RKJ,ZDKJ,KRKJ
C
      INTEGER IPIV(LINMAX*NATCLU_M),RHS_SLOT(LINMAX)
      INTEGER RHS_NRHS,IRHS,PACK_RHS,PACK_STATUS,L_RHS
      INTEGER NACTIVE,NTYP_COUNT
      INTEGER NCOMP,ICOMP
      INTEGER IBATCH,NBATCH,IBLOCK_SIZE,IBLOCK_START
C
      CHARACTER*1 CH
      CHARACTER*8 PACK_ENV
C
      DATA PI /3.1415926535898D0/
C
      ONEC=(1.D0,0.D0)
      IC=(0.D0,1.D0)
      ZEROC=(0.D0,0.D0)
      RHS_NRHS=0
      PACK_RHS=1
      PACK_ENV='1'
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_PACK_RHS',PACK_ENV,
     1                              STATUS=PACK_STATUS)
      IF(PACK_STATUS.EQ.0) THEN
        IF(PACK_ENV(1:1).EQ.'0') PACK_RHS=0
      ENDIF
      DO IRHS=1,LINMAX
        RHS_SLOT(IRHS)=0
      ENDDO
      IBESS=3
      CH='N'
C     The dense operator is square only at the active energy-dependent
C     angular cutoff.  Count its true dimension before allocating it;
C     LINMAX*NATCLU_M is only an upper bound and may be much larger.
      NACTIVE=0
      DO JTYP=1,N_PROT
        LMJ=LMAX(JTYP,JE)
        NTYP_COUNT=NATYP(JTYP)
        NACTIVE=NACTIVE+NTYP_COUNT*(LMJ+1)*(LMJ+1)
      ENDDO
C
C  Construction of the multiple scattering matrix MS = (I-GoT). 
C    Elements are stored using a linear index LINJ representing 
C    (J,LJ)
C
      ALLOCATE(SM(NACTIVE,NACTIVE))
      IF(NCOMP.GT.0) IN(1:NACTIVE,1:NCOMP)=ZEROC
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
            ILJ=LJ*LJ+LJ+1
            TLJ=DCMPLX(TL(LJ,1,JTYP,JE))
            DO MJ=-LJ,LJ
              INDJ=ILJ+MJ
              JLIN=JLIN+1
C
C             Pack the absorber channels selected by LF1:LF2:ISTEP_LF.
C             This skips angular-momentum channels that cannot contribute.
C
              IF((JTYP.EQ.1).AND.(JNUM.EQ.1)) THEN
                IF(PACK_RHS.EQ.1) THEN
                 IF((LJ.GE.LF1).AND.(LJ.LE.LF2)) THEN
                  IF(MOD(LJ-LF1,ISTEP_LF).EQ.0) THEN
                  RHS_NRHS=RHS_NRHS+1
                  RHS_SLOT(INDJ)=RHS_NRHS
                  ENDIF
                 ENDIF
                ELSEIF(LJ.LE.LF2) THEN
                  RHS_SLOT(INDJ)=JLIN
                ENDIF
              ENDIF
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
                    ATTKJ=DEXP(-DIMAG(DCMPLX(VK(JE)))*RKJ)
                    EXPKJ=(XKJ+IC*YKJ)/RKJ
                    ZDKJ=ZKJ/RKJ
                    CALL SPH_HAR2(2*NL_M,ZDKJ,EXPKJ,YLM,LMJ+LMK)
                    CALL BESPHE2(LMJ+LMK+1,IBESS,KRKJ,HL1)
                  ENDIF
C
                  DO LK=0,LMK
                    ILK=LK*LK+LK+1
                    L_MIN=ABS(LK-LJ)
                    L_MAX=LK+LJ
                    TLK=DCMPLX(TL(LK,1,KTYP,JE))
                    DO MK=-LK,LK
                      INDK=ILK+MK
                      KLIN=KLIN+1
                      SM(KLIN,JLIN)=ZEROC
                      SUM_L=ZEROC
                      IF(KATL.NE.JATL) THEN
                        CALL GAUNT2(LK,MK,LJ,MJ,GNT)
C
                        DO L=L_MIN,L_MAX,2
                          M=MJ-MK
                          IF(ABS(M).LE.L) THEN
                            SUM_L=SUM_L+(IC**L)*HL1(L)*
     1                                   YLM(L,M)*GNT(L)
                          ENDIF
                        ENDDO
                        SUM_L=SUM_L*ATTKJ*4.D0*PI*IC
                      ELSE
                        SUM_L=ZEROC
                      ENDIF
C
                      IF(KLIN.EQ.JLIN) THEN
                         SM(KLIN,JLIN)=ONEC-TLK*SUM_L
                         IF(NCOMP.GT.0) THEN
                           IF((JTYP.EQ.1).AND.(JNUM.EQ.1)) THEN
                             DO ICOMP=1,NCOMP
                               IN(KLIN,ICOMP)=SRC_COMP(INDJ,ICOMP)
                             ENDDO
                           ENDIF
                         ELSEIF((JTYP.EQ.1).AND.(JNUM.EQ.1).AND.
     1                          (RHS_SLOT(INDJ).GT.0)) THEN
                           IN(KLIN,RHS_SLOT(INDJ))=ONEC
                         ENDIF
                      ELSE
                         SM(KLIN,JLIN)=-TLK*SUM_L
                         IF(NCOMP.EQ.0.AND.(JTYP.EQ.1).AND.
     1                      (JNUM.EQ.1).AND.
     2                      (RHS_SLOT(INDJ).GT.0)) THEN
                           IN(KLIN,RHS_SLOT(INDJ))=ZEROC
                         ENDIF
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
C     Solve only the packed absorber channels selected by the PED final-L
C     range. Channels excluded by ISTEP_LF cannot contribute to the output.
      IF(PACK_RHS.EQ.1) THEN
        LW2=RHS_NRHS
      ELSE
        L_RHS=MIN(LF2,LMAX(1,JE))
        IF(L_RHS.LT.0) L_RHS=0
        LW2=(L_RHS+1)*(L_RHS+1)
      ENDIF
      IF(NCOMP.GT.0) LW2=NCOMP
      IF(NCOMP.GT.0) WRITE(6,*)
     1   'MSSPEC_SOURCE_COMPRESSION raw_rhs,combined_rhs=',
     2   RHS_NRHS,NCOMP
C
C   Partial inversion of the multiple scattering matrix MS and 
C     multiplication by T : the LAPACK subroutine performing 
C
C                          A * x = b
C
C     is used where b is the block column corresponding to 
C     the absorber 0 in the identity matrix. x is then TAU^{j 0}.
C
      CALL ZGETRF(JLIN,JLIN,SM,NACTIVE,IPIV,INFO1)
      IF(INFO1.NE.0) THEN
        WRITE(6,*) '     --->  INFO1 =',INFO1
      ELSE
        CALL ZGETRS(CH,JLIN,LW2,SM,NACTIVE,IPIV,
     1              IN,LINMAX*NATCLU_M,INFO)
      ENDIF
C
C   Storage of the Tau matrix
C
      JLIN=0
      DO JTYP=1,N_PROT
        NBTYPJ=NATYP(JTYP)
        LMJ=LMAX(JTYP,JE) 
        DO JNUM=1,NBTYPJ
          JATL=NCORR(JNUM,JTYP)
C
          DO LJ=0,LMJ
            ILJ=LJ*LJ+LJ+1
            TLJ=DCMPLX(TL(LJ,1,JTYP,JE))
            DO MJ=-LJ,LJ
              INDJ=ILJ+MJ
              JLIN=JLIN+1
C
              KLIN=0
              DO KTYP=1,N_PROT
                NBTYPK=NATYP(KTYP)
                LMK=LMAX(KTYP,JE)
                DO KNUM=1,NBTYPK
                  KATL=NCORR(KNUM,KTYP)
C
                  DO LK=0,LMK
                    ILK=LK*LK+LK+1
                    DO MK=-LK,LK
                      INDK=ILK+MK
                      KLIN=KLIN+1
                      IF(NCOMP.GT.0) THEN
                        IF((JATL.EQ.1).AND.(INDJ.EQ.1)) THEN
                          DO ICOMP=1,NCOMP
                            TAU(INDK,ICOMP,KATL)=
     1                           CMPLX(IN(KLIN,ICOMP))
                          ENDDO
                        ENDIF
                      ELSEIF((JATL.EQ.1).AND.
     1                       (RHS_SLOT(INDJ).GT.0)) THEN
                        TAU(INDK,INDJ,KATL)=
     1                    CMPLX(IN(KLIN,RHS_SLOT(INDJ))*TLJ)
                      ENDIF
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
      IF(ALLOCATED(SM)) DEALLOCATE(SM)
      RETURN
C
      END
