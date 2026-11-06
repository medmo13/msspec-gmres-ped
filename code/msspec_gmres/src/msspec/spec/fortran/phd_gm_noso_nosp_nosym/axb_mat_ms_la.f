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
      COMPLEX*16 XTEST(LINMAX*NATCLU_M)
      COMPLEX*16 YTEST(LINMAX*NATCLU_M),YREF(LINMAX*NATCLU_M)
      COMPLEX*16, ALLOCATABLE :: XBLOCK(:,:),BBLOCK(:,:)
      REAL*8, ALLOCATABLE :: RBLOCK(:)
      INTEGER, ALLOCATABLE :: ITBLOCK(:),INFOBLOCK(:)
C
      COMPLEX TAU(LINMAX,LINFMAX,NATCLU_M)
C
      REAL*8 PI,ATTKJ,GNT(0:N_GAUNT),XKJ,YKJ,ZKJ,RKJ,ZDKJ,KRKJ
C
      INTEGER IPIV(LINMAX*NATCLU_M),RHS_SLOT(LINMAX)
      INTEGER RHS_NRHS,IRHS,PACK_RHS,PACK_STATUS,L_RHS,CHECK_MV
      INTEGER NACTIVE,NTYP_COUNT
      INTEGER NCOMP,ICOMP
      INTEGER IBATCH,NBATCH,IBLOCK_SIZE,IBLOCK_START
      INTEGER USE_GMRES
      INTEGER USE_TREE_SERIES,TREE_ORDER,TREE_INFO
      REAL*8 ERR_MV,NRM_MV,MAX_SM,MAX_X,MAX_REF,MAX_TEST
C
      CHARACTER*1 CH
      CHARACTER*8 PACK_ENV
      CHARACTER*32 TREE_ENV
C
      DATA PI /3.1415926535898D0/
C
      ONEC=(1.D0,0.D0)
      IC=(0.D0,1.D0)
      ZEROC=(0.D0,0.D0)
      USE_GMRES=0
      USE_TREE_SERIES=0
      TREE_ENV='0'
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TEST_TREE_SERIES',
     1                              TREE_ENV,STATUS=TREE_INFO)
      IF(TREE_INFO.EQ.0) THEN
        IF(TREE_ENV(1:1).EQ.'1') THEN
          USE_TREE_SERIES=1
          USE_GMRES=1
        ENDIF
      ENDIF
      CHECK_MV=0
      PACK_ENV='0'
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_CHECK_MV',PACK_ENV,
     1                              STATUS=PACK_STATUS)
      IF(PACK_STATUS.EQ.0) THEN
        IF(PACK_ENV(1:1).EQ.'1') CHECK_MV=1
      ENDIF
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
      IF((USE_GMRES.EQ.0).OR.(CHECK_MV.EQ.1)) THEN
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
                      IF((USE_GMRES.EQ.0).OR.(CHECK_MV.EQ.1))
     1                  SM(KLIN,JLIN)=ZEROC
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
                         IF((USE_GMRES.EQ.0).OR.(CHECK_MV.EQ.1))
     1                     SM(KLIN,JLIN)=ONEC-TLK*SUM_L
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
                         IF((USE_GMRES.EQ.0).OR.(CHECK_MV.EQ.1))
     1                     SM(KLIN,JLIN)=-TLK*SUM_L
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
      ELSE
C       In matrix-free mode only form the absorber source vectors. Skip
C       evaluating every dense matrix element during setup.
        JLIN=0
        DO JTYP=1,N_PROT
          LMJ=LMAX(JTYP,JE)
          DO JNUM=1,NATYP(JTYP)
            DO LJ=0,LMJ
              ILJ=LJ*LJ+LJ+1
              DO MJ=-LJ,LJ
                JLIN=JLIN+1
                INDJ=ILJ+MJ
                IF((JTYP.EQ.1).AND.(JNUM.EQ.1)) THEN
                  IF(PACK_RHS.EQ.1) THEN
                    IF((LJ.GE.LF1).AND.(LJ.LE.LF2).AND.
     1                 (MOD(LJ-LF1,ISTEP_LF).EQ.0)) THEN
                      RHS_NRHS=RHS_NRHS+1
                      RHS_SLOT(INDJ)=RHS_NRHS
                    ENDIF
                  ELSEIF(LJ.LE.LF2) THEN
                    RHS_SLOT(INDJ)=JLIN
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IF(PACK_RHS.EQ.1) THEN
          LW2=RHS_NRHS
        ELSE
          L_RHS=MIN(LF2,LMAX(1,JE))
          IF(L_RHS.LT.0) L_RHS=0
          LW2=(L_RHS+1)*(L_RHS+1)
        ENDIF
        IF(NCOMP.GT.0) LW2=NCOMP
        DO IRHS=1,LW2
          IN(1:JLIN,IRHS)=ZEROC
        ENDDO
        JLIN=0
        DO JTYP=1,N_PROT
          LMJ=LMAX(JTYP,JE)
          DO JNUM=1,NATYP(JTYP)
            DO LJ=0,LMJ
              ILJ=LJ*LJ+LJ+1
              DO MJ=-LJ,LJ
                JLIN=JLIN+1
                INDJ=ILJ+MJ
                IF(NCOMP.GT.0) THEN
                  IF((JTYP.EQ.1).AND.(JNUM.EQ.1)) THEN
                    DO ICOMP=1,NCOMP
                      IN(JLIN,ICOMP)=SRC_COMP(INDJ,ICOMP)
                    ENDDO
                  ENDIF
                ELSEIF((JTYP.EQ.1).AND.(JNUM.EQ.1).AND.
     1                 (RHS_SLOT(INDJ).GT.0)) THEN
                  IN(JLIN,RHS_SLOT(INDJ))=ONEC
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
C     Optional exact-action check. It compares the streamed operator
C     against the dense matrix assembled above on one reproducible vector.
      IF(CHECK_MV.EQ.1) THEN
        DO I=1,JLIN
          XTEST(I)=DCMPLX(DSIN(DBLE(I)*0.173D0),
     1                    DCOS(DBLE(I)*0.119D0))
        ENDDO
        CALL ZGEMV('N',JLIN,JLIN,ONEC,SM,NACTIVE,
     1             XTEST,1,ZEROC,YREF,1)
        CALL MS_MATRIX_ACTION(JE,JLIN,XTEST,YTEST)
        ERR_MV=0.D0
        NRM_MV=0.D0
        DO I=1,JLIN
          ERR_MV=ERR_MV+ABS(YTEST(I)-YREF(I))**2
          NRM_MV=NRM_MV+ABS(YREF(I))**2
        ENDDO
        MAX_SM=MAXVAL(ABS(SM(1:JLIN,1:JLIN)))
        MAX_X=MAXVAL(ABS(XTEST(1:JLIN)))
        MAX_REF=MAXVAL(ABS(YREF(1:JLIN)))
        MAX_TEST=MAXVAL(ABS(YTEST(1:JLIN)))
        ERR_MV=DSQRT(ERR_MV)/MAX(DSQRT(NRM_MV),1.D-300)
        WRITE(6,*) 'MSSPEC_CHECK_MV relative error =',ERR_MV
        WRITE(6,*) 'MSSPEC_CHECK_MV max(SM,x,ref,stream)=',
     1             MAX_SM,MAX_X,MAX_REF,MAX_TEST
C       Check the cached block action against that same dense reference.
        ALLOCATE(XBLOCK(JLIN,1),BBLOCK(JLIN,1))
        XBLOCK(:,1)=XTEST(1:JLIN)
        CALL MS_MATRIX_ACTION_BLOCK(JE,JLIN,1,XBLOCK,BBLOCK)
        ERR_MV=DSQRT(SUM(ABS(BBLOCK(:,1)-YREF(1:JLIN))**2))
     1         /MAX(DSQRT(NRM_MV),1.D-300)
        WRITE(6,*) 'MSSPEC_CHECK_MV_BLOCK relative error =',ERR_MV
        DEALLOCATE(XBLOCK,BBLOCK)
      ENDIF
C
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
      GMRES_ALLOW_INEXACT=0
      IF(USE_TREE_SERIES.EQ.1) THEN
        TREE_ORDER=2
        TREE_ENV=' '
        CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_SERIES_ORDER',
     1                                TREE_ENV,STATUS=TREE_INFO)
        IF(TREE_INFO.EQ.0) THEN
          READ(TREE_ENV,*,IOSTAT=TREE_INFO) TREE_ORDER
          IF(TREE_INFO.NE.0) TREE_ORDER=2
        ENDIF
        IF(TREE_ORDER.LT.1.OR.TREE_ORDER.GT.200) THEN
          WRITE(6,*) 'MSSPEC_TREE_SERIES_ORDER must be 1..200'
          STOP 98
        ENDIF
        CALL MS_TREE_SERIES_BATCH(JE,JLIN,LW2,IN,TREE_ORDER,
     1                            TREE_INFO)
        IF(TREE_INFO.NE.0) THEN
          WRITE(6,*) 'MSSPEC_TEST_TREE_SERIES failed, info=',TREE_INFO
C         No dense fallback: a failed tree solve is a method failure.
          WRITE(6,*) 'MSSPEC_SOLVER_STATUS je,tree_info,status=',
     1               JE,TREE_INFO,'none_available'
          STOP 99
        ENDIF
      ENDIF
      IF(USE_GMRES.EQ.0) THEN
        CALL ZGETRF(JLIN,JLIN,SM,NACTIVE,IPIV,INFO1)
        IF(INFO1.NE.0) THEN
          WRITE(6,*) '     --->  INFO1 =',INFO1
        ELSE
          CALL ZGETRS(CH,JLIN,LW2,SM,NACTIVE,IPIV,
     1                IN,LINMAX*NATCLU_M,INFO)
        ENDIF
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

C=======================================================================
C
      SUBROUTINE MS_TREE_SERIES_BATCH(JE,N,NRHS,B,ORDER,INFO)
C
C   Exact block-tree solve.  Two solvers share the tree:
C     - default: finite preconditioned correction series (ORDER terms),
C       x <- A_T**(-1) (b - (A-A_T)x), used for the series-order study;
C     - MSSPEC_TREE_GMRES=1: restarted GMRES preconditioned by the tree.
C   The tree matrix contains the production diagonal and exact pair blocks.
C
      USE DIM_MOD
      USE COOR_MOD
      USE INIT_L_MOD
      USE TRANS_MOD
      INTEGER JE,N,NRHS,ORDER,INFO
      INTEGER I,J,K,L,T,NA,JTYP,JNUM,SITE,POS,QI,QP,QC,NMAX
      INTEGER, ALLOCATABLE :: STYPE(:),SNUM(:),Q(:),START(:)
      INTEGER, ALLOCATABLE :: PARENT(:),VISIT(:),IPIV(:,:),CHILD(:)
      COMPLEX*16, ALLOCATABLE :: LU(:,:,:),EP(:,:,:),EC(:,:,:)
      COMPLEX*16, ALLOCATABLE :: X(:,:),BSRC(:,:),FULL(:,:),TREE(:,:)
      COMPLEX*16, ALLOCATABLE :: WORK(:,:),TMP(:,:)
      COMPLEX*16, ALLOCATABLE :: DBLOCK(:,:)
      REAL*8, ALLOCATABLE :: BEST(:)
      REAL*8 D2,DX,DY,DZ,TREE_FROB2,OMIT_FROB2,TOTAL_FROB2
      REAL*8 TREE_TOL
      INTEGER INFO2,NR,USE_TREE_DIAG,DIAG_STATUS
      INTEGER USE_TREE_GMRES,GMRES_STATUS,GMRES_RESTART
      INTEGER GMRES_RESTART_STATUS,TOL_STATUS,VERIFY_STATUS
      INTEGER VERIFY_RES
      CHARACTER*32 TREE_DIAG_ENV,GMRES_ENV,GMRES_RESTART_ENV
      CHARACTER*32 TREE_TOL_ENV,VERIFY_ENV
      COMPLEX*16 B(LINMAX*NATCLU_M,NRHS)
      INTEGER*8 TC0,TC1,TC2,TCR
C
      CALL SYSTEM_CLOCK(TC0,TCR)
      INFO=0
      NA=0
      USE_TREE_DIAG=0
      USE_TREE_GMRES=0
      GMRES_ENV='0'
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_GMRES',
     1       GMRES_ENV,STATUS=GMRES_STATUS)
      IF(GMRES_STATUS.EQ.0) THEN
        IF(GMRES_ENV(1:1).EQ.'1') USE_TREE_GMRES=1
      ENDIF
      GMRES_RESTART=6
      GMRES_RESTART_ENV='6'
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_GMRES_RESTART',
     1       GMRES_RESTART_ENV,STATUS=GMRES_RESTART_STATUS)
      IF(GMRES_RESTART_STATUS.EQ.0) THEN
        READ(GMRES_RESTART_ENV,*,IOSTAT=INFO2) GMRES_RESTART
        IF(INFO2.NE.0) GMRES_RESTART=6
      ENDIF
      GMRES_RESTART=MAX(2,MIN(12,GMRES_RESTART))
      TREE_TOL=1.D-2
      TREE_TOL_ENV='1.D-2'
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_TOL',
     1       TREE_TOL_ENV,STATUS=TOL_STATUS)
      IF(TOL_STATUS.EQ.0) THEN
        READ(TREE_TOL_ENV,*,IOSTAT=INFO2) TREE_TOL
        IF(INFO2.NE.0) TREE_TOL=1.D-2
      ENDIF
      TREE_TOL=MAX(1.D-12,MIN(1.D-1,TREE_TOL))
      TREE_DIAG_ENV='0'
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_DIAGNOSTICS',
     1                              TREE_DIAG_ENV,STATUS=DIAG_STATUS)
      IF(DIAG_STATUS.EQ.0) THEN
        IF(TREE_DIAG_ENV(1:1).EQ.'1') USE_TREE_DIAG=1
      ENDIF
      DO JTYP=1,N_PROT
        NA=NA+NATYP(JTYP)
      ENDDO
      IF(NA.LE.0.OR.N.LE.0.OR.NRHS.LE.0) THEN
        INFO=-1
        RETURN
      ENDIF
      ALLOCATE(STYPE(NA),SNUM(NA),Q(NA),START(NA))
      ALLOCATE(PARENT(NA),VISIT(NA),CHILD(NA))
      ALLOCATE(BEST(NA))
      SITE=0
      POS=1
      DO JTYP=1,N_PROT
        DO JNUM=1,NATYP(JTYP)
          SITE=SITE+1
          STYPE(SITE)=JTYP
          SNUM(SITE)=JNUM
          Q(SITE)=(LMAX(JTYP,JE)+1)*(LMAX(JTYP,JE)+1)
          START(SITE)=POS
          POS=POS+Q(SITE)
        ENDDO
      ENDDO
      IF(POS-1.NE.N.OR.MAXVAL(Q).GT.LINMAX) THEN
        INFO=-2
        GOTO 900
      ENDIF
C
C   Prim's algorithm selects the geometric spanning tree (minimum distance).
C
      PARENT=0
      VISIT=0
      BEST=HUGE(1.D0)
      PARENT(1)=1
      BEST(1)=0.D0
      DO I=1,NA
        SITE=0
        D2=HUGE(1.D0)
        DO J=1,NA
          IF(VISIT(J).EQ.0.AND.BEST(J).LT.D2) THEN
            SITE=J
            D2=BEST(J)
          ENDIF
        ENDDO
        IF(SITE.EQ.0) THEN
          INFO=-3
          GOTO 900
        ENDIF
        VISIT(SITE)=I
        DO J=1,NA
          IF(VISIT(J).EQ.0) THEN
            DX=SYM_AT(1,NCORR(SNUM(SITE),STYPE(SITE)))
     1        -SYM_AT(1,NCORR(SNUM(J),STYPE(J)))
            DY=SYM_AT(2,NCORR(SNUM(SITE),STYPE(SITE)))
     1        -SYM_AT(2,NCORR(SNUM(J),STYPE(J)))
            DZ=SYM_AT(3,NCORR(SNUM(SITE),STYPE(SITE)))
     1        -SYM_AT(3,NCORR(SNUM(J),STYPE(J)))
            D2=DX*DX+DY*DY+DZ*DZ
            IF(D2.LT.BEST(J)) THEN
              BEST(J)=D2
              PARENT(J)=SITE
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      DO I=1,NA
        CHILD(VISIT(I))=I
      ENDDO
      VISIT(1:NA)=CHILD(1:NA)
C
C   EP(:,:,child) is A(parent,child); EC is A(child,parent).
C
      ALLOCATE(EP(LINMAX,LINMAX,NA),EC(LINMAX,LINMAX,NA))
      ALLOCATE(LU(LINMAX,LINMAX,NA),IPIV(LINMAX,NA))
      EP=DCMPLX(0.D0,0.D0)
      EC=DCMPLX(0.D0,0.D0)
      LU=DCMPLX(0.D0,0.D0)
      DO I=1,NA
        QI=Q(I)
        DO J=1,QI
          LU(J,J,I)=DCMPLX(1.D0,0.D0)
        ENDDO
      ENDDO
      DO I=2,NA
        SITE=PARENT(I)
        CALL MS_TREE_PAIR_BLOCK(JE,STYPE(SITE),SNUM(SITE),
     1       STYPE(I),SNUM(I),EP(1,1,I))
        CALL MS_TREE_PAIR_BLOCK(JE,STYPE(I),SNUM(I),
     1       STYPE(SITE),SNUM(SITE),EC(1,1,I))
      ENDDO
C
C   Optional exact Frobenius norms of retained and omitted directed blocks.
C   This adds one full pair-block traversal and is for diagnostics only.
C
      IF(USE_TREE_DIAG.EQ.1) THEN
        TREE_FROB2=0.D0
        OMIT_FROB2=0.D0
        DO I=2,NA
          SITE=PARENT(I)
          QP=Q(SITE)
          QC=Q(I)
          TREE_FROB2=TREE_FROB2+
     1       SUM(ABS(EP(1:QP,1:QC,I))**2)+
     2       SUM(ABS(EC(1:QC,1:QP,I))**2)
        ENDDO
        ALLOCATE(DBLOCK(LINMAX,LINMAX))
        DO I=1,NA
          DO J=1,NA
            IF(I.NE.J.AND.PARENT(I).NE.J.AND.PARENT(J).NE.I) THEN
              CALL MS_TREE_PAIR_BLOCK(JE,STYPE(I),SNUM(I),
     1             STYPE(J),SNUM(J),DBLOCK)
              OMIT_FROB2=OMIT_FROB2+
     1             SUM(ABS(DBLOCK(1:Q(I),1:Q(J)))**2)
            ENDIF
          ENDDO
        ENDDO
        DEALLOCATE(DBLOCK)
        TOTAL_FROB2=TREE_FROB2+OMIT_FROB2
        WRITE(6,*) 'TREE_BLOCK_NORMS',DSQRT(TREE_FROB2),
     1       DSQRT(OMIT_FROB2),DSQRT(TOTAL_FROB2),
     2       DSQRT(TREE_FROB2/MAX(TOTAL_FROB2,1.D-300)),
     3       DSQRT(OMIT_FROB2/MAX(TOTAL_FROB2,1.D-300))
      ENDIF
C
C   Eliminate leaves toward the root, factoring each diagonal Schur block.
C
      ALLOCATE(TMP(LINMAX,LINMAX))
      DO K=NA,2,-1
        I=VISIT(K)
        QI=Q(I)
        DO J=2,NA
          IF(PARENT(J).EQ.I) THEN
            QC=Q(J)
            QP=QI
            TMP(1:QC,1:QP)=EC(1:QC,1:QP,J)
            CALL ZGETRS('N',QC,QP,LU(1,1,J),LINMAX,
     1                  IPIV(1,J),TMP,LINMAX,INFO2)
            IF(INFO2.NE.0) THEN
              INFO=-10-INFO2
              GOTO 900
            ENDIF
            DO T=1,QP
              DO POS=1,QP
              DO L=1,QC
                  LU(POS,T,I)=LU(POS,T,I)-
     1              EP(POS,L,J)*TMP(L,T)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        CALL ZGETRF(QI,QI,LU(1,1,I),LINMAX,IPIV(1,I),INFO2)
        IF(INFO2.NE.0) THEN
          INFO=-20-INFO2
          GOTO 900
        ENDIF
      ENDDO
      CALL ZGETRF(Q(1),Q(1),LU(1,1,1),LINMAX,IPIV(1,1),INFO2)
      IF(INFO2.NE.0) THEN
        INFO=-30-INFO2
        GOTO 900
      ENDIF
C
C   x0 = inverse(tree) b.
C
      CALL SYSTEM_CLOCK(TC1,TCR)
      NMAX=LINMAX*NATCLU_M
      ALLOCATE(BSRC(NMAX,NRHS),FULL(NMAX,NRHS))
      BSRC=DCMPLX(0.D0,0.D0)
      BSRC(1:N,1:NRHS)=B(1:N,1:NRHS)
      IF(USE_TREE_GMRES.EQ.1) THEN
        CALL MS_TREE_GMRES(JE,N,NA,NRHS,Q,START,PARENT,VISIT,
     1       LU,IPIV,EP,EC,BSRC,B,ORDER,TREE_TOL,GMRES_RESTART,INFO2)
        IF(INFO2.NE.0) THEN
          INFO=INFO2
          GOTO 900
        ENDIF
        GOTO 880
      ENDIF
      ALLOCATE(X(NMAX,NRHS),TREE(NMAX,NRHS),WORK(NMAX,NRHS))
      X=BSRC
      CALL MS_TREE_SOLVE(NA,NRHS,Q,START,PARENT,VISIT,
     1     LU,IPIV,EP,EC,X,INFO)
      IF(INFO.NE.0) GOTO 900
C
C   Each correction uses the exact production A action; only its tree part
C   is stored. This costs matrix-free A actions, not dense matrix storage.
C
      DO NR=1,ORDER
        FULL=DCMPLX(0.D0,0.D0)
        CALL MS_MATRIX_ACTION_BLOCK(JE,N,NRHS,X,FULL)
        TREE=X
        DO I=2,NA
          SITE=PARENT(I)
          QP=Q(SITE)
          QC=Q(I)
          DO T=1,NRHS
            DO J=1,QP
              DO K=1,QC
                TREE(START(SITE)+J-1,T)=
     1             TREE(START(SITE)+J-1,T)+
     2             EP(J,K,I)*X(START(I)+K-1,T)
              ENDDO
            ENDDO
            DO J=1,QC
              DO K=1,QP
                TREE(START(I)+J-1,T)=TREE(START(I)+J-1,T)+
     1             EC(J,K,I)*X(START(SITE)+K-1,T)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        WORK=BSRC-(FULL-TREE)
        CALL MS_TREE_SOLVE(NA,NRHS,Q,START,PARENT,VISIT,
     1       LU,IPIV,EP,EC,WORK,INFO)
        IF(INFO.NE.0) GOTO 900
        X=WORK
      ENDDO
      B(1:N,1:NRHS)=X(1:N,1:NRHS)
C
  880 CONTINUE
      CALL SYSTEM_CLOCK(TC2,TCR)
      WRITE(6,*) 'MSSPEC_TIMING je,tree_setup_s,solve_s,n,na=',JE,
     1     DBLE(TC1-TC0)/DBLE(TCR),DBLE(TC2-TC1)/DBLE(TCR),N,NA
C   The extra full-operator matvec is opt-in (MSSPEC_TREE_VERIFY=1 or
C   MSSPEC_TREE_DIAGNOSTICS=1); the GMRES stopping test already
C   certifies the preconditioned residual.
      VERIFY_RES=USE_TREE_DIAG
      VERIFY_ENV=' '
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_VERIFY',VERIFY_ENV,
     1                              STATUS=VERIFY_STATUS)
      IF(VERIFY_STATUS.EQ.0.AND.VERIFY_ENV(1:1).EQ.'1') VERIFY_RES=1
      IF(VERIFY_RES.EQ.1) THEN
      FULL=DCMPLX(0.D0,0.D0)
      CALL MS_MATRIX_ACTION_BLOCK(JE,N,NRHS,B,FULL)
      DO T=1,NRHS
        DX=DSQRT(SUM(ABS(FULL(1:N,T)-BSRC(1:N,T))**2))
        DY=DSQRT(SUM(ABS(BSRC(1:N,T))**2))
        WRITE(6,*) 'TREE_SERIES order,rhs,res=',ORDER,T,
     1             DX/MAX(DY,1.D-300)
      ENDDO
      ENDIF
  900 CONTINUE
      IF(ALLOCATED(X)) DEALLOCATE(X)
      IF(ALLOCATED(BSRC)) DEALLOCATE(BSRC)
      IF(ALLOCATED(FULL)) DEALLOCATE(FULL)
      IF(ALLOCATED(TREE)) DEALLOCATE(TREE)
      IF(ALLOCATED(WORK)) DEALLOCATE(WORK)
      IF(ALLOCATED(EP)) DEALLOCATE(EP,EC,LU,IPIV,TMP)
      IF(ALLOCATED(STYPE)) DEALLOCATE(STYPE,SNUM,Q,START,
     1                                 PARENT,VISIT,CHILD,BEST)
      RETURN
      END

C=======================================================================
C
      SUBROUTINE MS_TREE_GMRES(JE,N,NA,NRHS,Q,START,PARENT,VISIT,
     1       LU,IPIV,EP,EC,RHS,XOUT,MAXIT,TOL,RESTART,INFO)
C
C   Restarted GMRES on the left-preconditioned system M**(-1) A x =
C   M**(-1) b, where M is the existing exact block-tree factorization.
C   A is applied by the production streamed block action.
C
C   Each Arnoldi column is reduced at once by complex Givens rotations.
C   After column J the rotated right-hand side gives the preconditioned
C   residual norm ABS(G(J+1)) without forming the solution. A right-hand
C   side leaves the inner loop as soon as that estimate meets TOL (early
C   stop). The true preconditioned residual is recomputed at the end of
C   each cycle, and only that residual decides convergence.
C
      USE DIM_MOD
      USE, INTRINSIC :: IEEE_ARITHMETIC
      INTEGER JE,N,NA,NRHS,Q(NA),START(NA),PARENT(NA),VISIT(NA)
      INTEGER MAXIT,RESTART,INFO,IPIV(LINMAX,NA)
      INTEGER NMAX,MMAX,I,J,L,IRHS,NACT,IDX,INFO2
      INTEGER, ALLOCATABLE :: ACTIVE(:),INNER(:),K(:),ITERS(:),AIDX(:)
      INTEGER, ALLOCATABLE :: BROKE(:)
      COMPLEX*16 LU(LINMAX,LINMAX,NA),EP(LINMAX,LINMAX,NA)
      COMPLEX*16 EC(LINMAX,LINMAX,NA)
      COMPLEX*16 RHS(LINMAX*NATCLU_M,NRHS)
      COMPLEX*16 XOUT(LINMAX*NATCLU_M,NRHS)
      COMPLEX*16, ALLOCATABLE :: V(:,:,:),H(:,:,:),G(:,:),SN(:,:)
      COMPLEX*16, ALLOCATABLE :: BPRE(:,:),R(:,:),PACK(:,:)
      COMPLEX*16, ALLOCATABLE :: W(:,:),AX(:,:)
      COMPLEX*16 ZDOTC,HA,HB,TOP
      REAL*8 TOL,BETA,WNORM,RAD,RESID
      REAL*8, ALLOCATABLE :: BNORMS(:),RNORMS(:),CR(:,:)
      EXTERNAL ZDOTC
      INTEGER ITER_DUMP,DUMP_STATUS,IU,STOPK,MAXEFF
      CHARACTER*8 DUMP_ENV,STOP_ENV
      INFO=0
      ITER_DUMP=0
      IU=77
      DUMP_ENV='0'
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_ITER_DUMP',
     1                              DUMP_ENV,STATUS=DUMP_STATUS)
      IF(DUMP_STATUS.EQ.0) THEN
        IF(DUMP_ENV(1:1).EQ.'1') THEN
          OPEN(UNIT=IU,FILE='msspec_tree_iter.dat',STATUS='REPLACE',
     1         FORM='FORMATTED',ACTION='WRITE',IOSTAT=INFO2)
          IF(INFO2.EQ.0) ITER_DUMP=1
        ENDIF
      ENDIF
C     MSSPEC_TREE_STOP_AFTER=k (study only, E16): return the iterate after
C     exactly k Arnoldi steps instead of solving to TOL. k=0 returns the
C     preconditioner alone, M**(-1) b. Unset or -1: normal solve.
C     A value that is not an integer, or k above the restart length, is
C     an error (INFO=-210), never a silent normal solve.
      STOPK=-1
      STOP_ENV='-1'
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_STOP_AFTER',
     1                              STOP_ENV,STATUS=DUMP_STATUS)
      IF(DUMP_STATUS.EQ.0) THEN
        READ(STOP_ENV,*,IOSTAT=INFO2) STOPK
        IF(INFO2.NE.0) STOPK=-2
      ENDIF
      MAXEFF=MAXIT
      IF(STOPK.GT.0) MAXEFF=STOPK
      NMAX=LINMAX*NATCLU_M
      MMAX=MAX(2,MIN(12,RESTART))
      ALLOCATE(V(NMAX,MMAX+1,NRHS),H(MMAX+1,MMAX,NRHS))
      ALLOCATE(G(MMAX+1,NRHS),SN(MMAX,NRHS),CR(MMAX,NRHS))
      ALLOCATE(BPRE(NMAX,NRHS),R(NMAX,NRHS))
      ALLOCATE(PACK(NMAX,NRHS),W(NMAX,NRHS),AX(NMAX,NRHS))
      ALLOCATE(ACTIVE(NRHS),INNER(NRHS),K(NRHS),ITERS(NRHS))
      ALLOCATE(AIDX(NRHS),BROKE(NRHS),BNORMS(NRHS),RNORMS(NRHS))
      BPRE=RHS
      CALL MS_TREE_SOLVE(NA,NRHS,Q,START,PARENT,VISIT,
     1     LU,IPIV,EP,EC,BPRE,INFO2)
      IF(INFO2.NE.0) THEN
        INFO=INFO2
        GOTO 990
      ENDIF
      XOUT=DCMPLX(0.D0,0.D0)
      R=BPRE
      DO IRHS=1,NRHS
        BNORMS(IRHS)=DSQRT(SUM(ABS(BPRE(1:N,IRHS))**2))
        ITERS(IRHS)=0
        ACTIVE(IRHS)=1
        BROKE(IRHS)=0
        IF(.NOT.IEEE_IS_FINITE(BNORMS(IRHS)) ) THEN
          INFO=-150-IRHS
          GOTO 990
        ELSEIF(BNORMS(IRHS).LE.1.D-280) THEN
          ACTIVE(IRHS)=0
        ENDIF
      ENDDO
      IF(ITER_DUMP.EQ.1) THEN
C       Iterate 0 is the preconditioner alone, M**(-1) b (BPRE), written
C       before the first Arnoldi step. Its residual is the absolute
C       preconditioned one, the same unit as ABS(G(J+1)) in the later
C       ITERATE blocks. AX is scratch here and is overwritten by the
C       first cycle, so the solve is unchanged.
        CALL MS_MATRIX_ACTION_BLOCK(JE,N,NRHS,BPRE,AX)
        CALL MS_TREE_SOLVE(NA,NRHS,Q,START,PARENT,VISIT,
     1       LU,IPIV,EP,EC,AX,INFO2)
        IF(INFO2.NE.0) THEN
          INFO=INFO2
          GOTO 990
        ENDIF
        DO IRHS=1,NRHS
          RESID=DSQRT(SUM(ABS(BPRE(1:N,IRHS)-AX(1:N,IRHS))**2))
          WRITE(IU,'(A,I6,I8,I6,I8,ES16.8)') 'ITERATE',IRHS,0,0,N,RESID
          DO I=1,N
            WRITE(IU,'(1P,2E24.16)') DBLE(BPRE(I,IRHS)),
     1           DIMAG(BPRE(I,IRHS))
          ENDDO
        ENDDO
      ENDIF
      IF(STOPK.LT.-1.OR.STOPK.GT.MMAX) THEN
        INFO=-210
        GOTO 990
      ENDIF
      IF(STOPK.GE.0) WRITE(6,*) 'TREE_GMRES stop_after=',STOPK
      IF(STOPK.EQ.0) THEN
        XOUT(1:N,1:NRHS)=BPRE(1:N,1:NRHS)
        ACTIVE=0
      ENDIF
      DO WHILE(ANY(ACTIVE.EQ.1))
        H=DCMPLX(0.D0,0.D0)
        G=DCMPLX(0.D0,0.D0)
        V=DCMPLX(0.D0,0.D0)
        K=0
        INNER=0
        DO IRHS=1,NRHS
          IF(ACTIVE(IRHS).EQ.1) THEN
            BETA=DSQRT(SUM(ABS(R(1:N,IRHS))**2))
            IF(.NOT.IEEE_IS_FINITE(BETA)) THEN
              INFO=-160-IRHS
              GOTO 990
            ELSEIF(BETA.LE.TOL*MAX(BNORMS(IRHS),1.D-280)) THEN
              ACTIVE(IRHS)=0
            ELSE
              V(1:N,1,IRHS)=R(1:N,IRHS)/BETA
              G(1,IRHS)=DCMPLX(BETA,0.D0)
              INNER(IRHS)=1
            ENDIF
          ENDIF
        ENDDO
        IF(ALL(ACTIVE.EQ.0)) EXIT
        DO J=1,MMAX
          NACT=0
          DO IRHS=1,NRHS
            IF(INNER(IRHS).EQ.1.AND.ITERS(IRHS).LT.MAXEFF) THEN
              NACT=NACT+1
              AIDX(NACT)=IRHS
              PACK(1:N,NACT)=V(1:N,J,IRHS)
            ENDIF
          ENDDO
          IF(NACT.EQ.0) EXIT
          CALL MS_MATRIX_ACTION_BLOCK(JE,N,NACT,PACK,W)
          CALL MS_TREE_SOLVE(NA,NACT,Q,START,PARENT,VISIT,
     1         LU,IPIV,EP,EC,W,INFO2)
          IF(INFO2.NE.0) THEN
            INFO=INFO2
            GOTO 990
          ENDIF
          DO IDX=1,NACT
            IRHS=AIDX(IDX)
            DO I=1,J
              H(I,J,IRHS)=ZDOTC(N,V(1,I,IRHS),1,W(1,IDX),1)
              W(1:N,IDX)=W(1:N,IDX)-H(I,J,IRHS)*
     1                             V(1:N,I,IRHS)
            ENDDO
            WNORM=DSQRT(SUM(ABS(W(1:N,IDX))**2))
            IF(.NOT.IEEE_IS_FINITE(WNORM)) THEN
              INFO=-170-IRHS
              GOTO 990
            ENDIF
            H(J+1,J,IRHS)=DCMPLX(WNORM,0.D0)
C           Earlier rotations act on the new column first.
            DO I=1,J-1
              HA=CR(I,IRHS)*H(I,J,IRHS)+SN(I,IRHS)*H(I+1,J,IRHS)
              HB=-DCONJG(SN(I,IRHS))*H(I,J,IRHS)+
     1            CR(I,IRHS)*H(I+1,J,IRHS)
              H(I,J,IRHS)=HA
              H(I+1,J,IRHS)=HB
            ENDDO
C           New rotation annihilates H(J+1,J).
            HA=H(J,J,IRHS)
            HB=H(J+1,J,IRHS)
            RAD=DSQRT(ABS(HA)**2+ABS(HB)**2)
            IF(ABS(HA).EQ.0.D0) THEN
              CR(J,IRHS)=0.D0
              SN(J,IRHS)=DCMPLX(1.D0,0.D0)
              TOP=HB
            ELSE
              CR(J,IRHS)=ABS(HA)/RAD
              SN(J,IRHS)=(HA/ABS(HA))*DCONJG(HB)/RAD
              TOP=(HA/ABS(HA))*RAD
            ENDIF
            H(J,J,IRHS)=TOP
            H(J+1,J,IRHS)=DCMPLX(0.D0,0.D0)
            HA=CR(J,IRHS)*G(J,IRHS)+SN(J,IRHS)*G(J+1,IRHS)
            G(J+1,IRHS)=-DCONJG(SN(J,IRHS))*G(J,IRHS)+
     1                  CR(J,IRHS)*G(J+1,IRHS)
            G(J,IRHS)=HA
            K(IRHS)=J
            ITERS(IRHS)=ITERS(IRHS)+1
            IF(STOPK.LT.0.AND.ABS(G(J+1,IRHS)).LE.
     1         TOL*MAX(BNORMS(IRHS),1.D-280)) THEN
              INNER(IRHS)=0
            ELSEIF(WNORM.GT.1.D-14) THEN
              V(1:N,J+1,IRHS)=W(1:N,IDX)/WNORM
            ELSE
              INNER(IRHS)=0
              BROKE(IRHS)=1
            ENDIF
            IF(ITER_DUMP.EQ.1) CALL MS_TREE_ITER_DUMP(N,J,IRHS,
     1           ITERS(IRHS),ABS(G(J+1,IRHS)),XOUT,V,H,G,
     1           NMAX,MMAX,NRHS,IU)
          ENDDO
        ENDDO
        DO IRHS=1,NRHS
          IF(K(IRHS).GT.0) THEN
C           Back substitution on the rotated triangular system.
            DO I=K(IRHS),1,-1
              HA=G(I,IRHS)
              DO L=I+1,K(IRHS)
                HA=HA-H(I,L,IRHS)*G(L,IRHS)
              ENDDO
              IF(ABS(H(I,I,IRHS)).EQ.0.D0) THEN
                INFO=-180-IRHS
                GOTO 990
              ENDIF
              G(I,IRHS)=HA/H(I,I,IRHS)
            ENDDO
            DO J=1,K(IRHS)
              XOUT(1:N,IRHS)=XOUT(1:N,IRHS)+
     1             V(1:N,J,IRHS)*G(J,IRHS)
            ENDDO
          ENDIF
        ENDDO
C       Stop mode: the iterate after STOPK steps is the result. No
C       residual check and no convergence error.
        IF(STOPK.GT.0) THEN
          ACTIVE=0
          BROKE=0
          EXIT
        ENDIF
        CALL MS_MATRIX_ACTION_BLOCK(JE,N,NRHS,XOUT,AX)
        CALL MS_TREE_SOLVE(NA,NRHS,Q,START,PARENT,VISIT,
     1       LU,IPIV,EP,EC,AX,INFO2)
        IF(INFO2.NE.0) THEN
          INFO=INFO2
          GOTO 990
        ENDIF
        R=BPRE-AX
        DO IRHS=1,NRHS
          IF(ACTIVE(IRHS).EQ.1.OR.BROKE(IRHS).EQ.1) THEN
            RNORMS(IRHS)=DSQRT(SUM(ABS(R(1:N,IRHS))**2))
            IF(.NOT.IEEE_IS_FINITE(RNORMS(IRHS))) THEN
              INFO=-190-IRHS
              GOTO 990
            ELSEIF(RNORMS(IRHS).LE.
     1             TOL*MAX(BNORMS(IRHS),1.D-280)) THEN
              ACTIVE(IRHS)=0
              BROKE(IRHS)=0
            ELSEIF(BROKE(IRHS).EQ.1.OR.ITERS(IRHS).GE.MAXIT) THEN
              INFO=-200-IRHS
              GOTO 990
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      DO IRHS=1,NRHS
        WRITE(6,*) 'TREE_GMRES iterations,rhs=',ITERS(IRHS),IRHS
      ENDDO
C     Final solution, written after the last cycle so the test can check
C     that the last ITERATE block of each RHS equals the returned X.
      IF(ITER_DUMP.EQ.1) THEN
        DO IRHS=1,NRHS
          WRITE(IU,'(A,I6,I8)') 'FINAL',IRHS,N
          DO I=1,N
            WRITE(IU,'(1P,2E24.16)') DBLE(XOUT(I,IRHS)),
     1           DIMAG(XOUT(I,IRHS))
          ENDDO
        ENDDO
      ENDIF
  990 CONTINUE
      IF(ITER_DUMP.EQ.1) CLOSE(IU)
      DEALLOCATE(V,H,G,SN,CR,BPRE,R,PACK,W,AX,ACTIVE,INNER,K,ITERS,
     1           AIDX,BROKE,BNORMS,RNORMS)
      RETURN
      END


C=======================================================================
C
C     Writes the partial GMRES iterate x_j = x0 + V_j y_j of one
C     right-hand side after Arnoldi step J. The triangular system is
C     solved into a local copy, so the solver state is not changed.
C     Runs only when MSSPEC_TREE_ITER_DUMP=1; output goes to unit IU.
      SUBROUTINE MS_TREE_ITER_DUMP(N,J,IRHS,ITER,RESID,XBASE,
     1           V,H,G,NMAX,MMAX,NRHS,IU)
      INTEGER N,J,IRHS,ITER,NMAX,MMAX,NRHS,IU,I,L
      REAL*8 RESID
      COMPLEX*16 XBASE(NMAX,NRHS)
      COMPLEX*16 V(NMAX,MMAX+1,NRHS),H(MMAX+1,MMAX,NRHS)
      COMPLEX*16 G(MMAX+1,NRHS)
      COMPLEX*16 YT(MMAX),XT(NMAX),HA
      DO I=J,1,-1
        HA=G(I,IRHS)
        DO L=I+1,J
          HA=HA-H(I,L,IRHS)*YT(L)
        ENDDO
        IF(ABS(H(I,I,IRHS)).EQ.0.D0) RETURN
        YT(I)=HA/H(I,I,IRHS)
      ENDDO
      XT(1:N)=XBASE(1:N,IRHS)
      DO L=1,J
        XT(1:N)=XT(1:N)+V(1:N,L,IRHS)*YT(L)
      ENDDO
      WRITE(IU,'(A,I6,I8,I6,I8,ES16.8)') 'ITERATE',IRHS,ITER,J,N,RESID
      DO I=1,N
        WRITE(IU,'(1P,2E24.16)') DBLE(XT(I)),DIMAG(XT(I))
      ENDDO
      RETURN
      END


C=======================================================================
C
      SUBROUTINE MS_TREE_SOLVE(NA,NRHS,Q,START,PARENT,VISIT,
     1          LU,IPIV,EP,EC,B,INFO)
      USE DIM_MOD
      INTEGER NA,NRHS,Q(NA),START(NA),PARENT(NA),VISIT(NA)
      INTEGER IPIV(LINMAX,NA),INFO,I,J,K,L,T,QP,QC,INFO2
      COMPLEX*16 LU(LINMAX,LINMAX,NA),EP(LINMAX,LINMAX,NA)
      COMPLEX*16 EC(LINMAX,LINMAX,NA)
      COMPLEX*16 B(LINMAX*NATCLU_M,NRHS)
      COMPLEX*16, ALLOCATABLE :: TMP(:,:)
      LOGICAL, SAVE :: NOPRE_INIT=.FALSE.,NOPRE=.FALSE.
      CHARACTER*8 NOPRE_ENV
      INTEGER NOPRE_ST
      INFO=0
C     MSSPEC_TREE_NOPRECOND=1 turns the tree solve into the identity
C     (M=I, GMRES without preconditioner); used only for the ablation.
      IF(.NOT.NOPRE_INIT) THEN
        NOPRE_INIT=.TRUE.
        NOPRE_ENV='0'
        CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_NOPRECOND',
     1       NOPRE_ENV,STATUS=NOPRE_ST)
        IF(NOPRE_ST.EQ.0.AND.NOPRE_ENV(1:1).EQ.'1') NOPRE=.TRUE.
      ENDIF
      IF(NOPRE) RETURN
      ALLOCATE(TMP(LINMAX,NRHS))
C     Child elimination: update parent RHS after solving child block.
      DO K=NA,2,-1
        I=VISIT(K)
        QP=Q(PARENT(I))
        QC=Q(I)
        TMP=DCMPLX(0.D0,0.D0)
        TMP(1:QC,1:NRHS)=B(START(I):START(I)+QC-1,1:NRHS)
        CALL ZGETRS('N',QC,NRHS,LU(1,1,I),LINMAX,
     1              IPIV(1,I),TMP,LINMAX,INFO2)
        IF(INFO2.NE.0) THEN
          INFO=-40-INFO2
          DEALLOCATE(TMP)
          RETURN
        ENDIF
        DO T=1,NRHS
          DO J=1,QP
            DO L=1,QC
              B(START(PARENT(I))+J-1,T)=
     1          B(START(PARENT(I))+J-1,T)-EP(J,L,I)*TMP(L,T)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C     Root solve followed by forward substitution through the tree.
      TMP=DCMPLX(0.D0,0.D0)
      TMP(1:Q(1),1:NRHS)=B(START(1):START(1)+Q(1)-1,1:NRHS)
      CALL ZGETRS('N',Q(1),NRHS,LU(1,1,1),LINMAX,
     1            IPIV(1,1),TMP,LINMAX,INFO2)
      IF(INFO2.NE.0) INFO=-50-INFO2
      B(START(1):START(1)+Q(1)-1,1:NRHS)=TMP(1:Q(1),1:NRHS)
      DO K=2,NA
        I=VISIT(K)
        QP=Q(PARENT(I))
        QC=Q(I)
        DO T=1,NRHS
          DO J=1,QC
            DO L=1,QP
              B(START(I)+J-1,T)=B(START(I)+J-1,T)-
     1          EC(J,L,I)*B(START(PARENT(I))+L-1,T)
            ENDDO
          ENDDO
        ENDDO
        TMP=DCMPLX(0.D0,0.D0)
        TMP(1:QC,1:NRHS)=B(START(I):START(I)+QC-1,1:NRHS)
        CALL ZGETRS('N',QC,NRHS,LU(1,1,I),LINMAX,
     1              IPIV(1,I),TMP,LINMAX,INFO2)
        IF(INFO2.NE.0) INFO=-60-INFO2
        B(START(I):START(I)+QC-1,1:NRHS)=TMP(1:QC,1:NRHS)
      ENDDO
      DEALLOCATE(TMP)
      RETURN
      END

C=======================================================================
C
      SUBROUTINE MS_TREE_PAIR_BLOCK(JE,DTYPE,DNUM,STYPE,SNUM,B)
      USE DIM_MOD
      USE COOR_MOD
      USE INIT_L_MOD
      USE TRANS_MOD
      INTEGER JE,DTYPE,DNUM,STYPE,SNUM,LJ,LK,MJ,MK,LD,LS
      INTEGER LMIN,LMAXV,L,M,IBESS,ID,IS
      COMPLEX*16 B(LINMAX,LINMAX),HL1(0:2*NL_M)
      COMPLEX*16 SUM_L,IC,EXPKJ,TLK,YLM(0:2*NL_M,-2*NL_M:2*NL_M)
      REAL*8 PI,ATTKJ,XKJ,YKJ,ZKJ,RKJ,ZDKJ,KRKJ,GNT(0:N_GAUNT)
      DATA PI /3.1415926535898D0/
      IC=DCMPLX(0.D0,1.D0)
      B=DCMPLX(0.D0,0.D0)
      XKJ=DBLE(SYM_AT(1,NCORR(DNUM,DTYPE))-
     1         SYM_AT(1,NCORR(SNUM,STYPE)))
      YKJ=DBLE(SYM_AT(2,NCORR(DNUM,DTYPE))-
     1         SYM_AT(2,NCORR(SNUM,STYPE)))
      ZKJ=DBLE(SYM_AT(3,NCORR(DNUM,DTYPE))-
     1         SYM_AT(3,NCORR(SNUM,STYPE)))
      RKJ=DSQRT(XKJ*XKJ+YKJ*YKJ+ZKJ*ZKJ)
      KRKJ=DBLE(VK(JE))*RKJ
      ATTKJ=DEXP(-DIMAG(DCMPLX(VK(JE)))*RKJ)
      EXPKJ=(XKJ+IC*YKJ)/RKJ
      ZDKJ=ZKJ/RKJ
      LD=LMAX(DTYPE,JE)
      LS=LMAX(STYPE,JE)
      CALL SPH_HAR2(2*NL_M,ZDKJ,EXPKJ,YLM,LD+LS)
      CALL BESPHE2(LD+LS+1,3,KRKJ,HL1)
      IS=0
      DO LJ=0,LS
        DO MJ=-LJ,LJ
          IS=IS+1
          ID=0
          DO LK=0,LD
            LMIN=ABS(LK-LJ)
            LMAXV=LK+LJ
            TLK=DCMPLX(TL(LK,1,DTYPE,JE))
            DO MK=-LK,LK
              ID=ID+1
              CALL GAUNT2(LK,MK,LJ,MJ,GNT)
              SUM_L=DCMPLX(0.D0,0.D0)
              DO L=LMIN,LMAXV,2
                M=MJ-MK
                IF(ABS(M).LE.L) SUM_L=SUM_L+(IC**L)*HL1(L)*
     1                                      YLM(L,M)*GNT(L)
              ENDDO
              SUM_L=SUM_L*ATTKJ*4.D0*PI*IC
              B(ID,IS)=-TLK*SUM_L
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END

C=======================================================================
C
      SUBROUTINE MS_MATRIX_ACTION(JE,N,X,Y)
C
C   Stream the production multiple-scattering matrix action Y=(I-GoT)X.
C   The element expression matches INV_MAT_MS. No dense matrix is stored.
C
      USE DIM_MOD
      USE COOR_MOD
      USE INIT_L_MOD
      USE TRANS_MOD
C
      COMPLEX*16 X(LINMAX*NATCLU_M),Y(LINMAX*NATCLU_M)
      COMPLEX*16 HL1(0:2*NL_M),SUM_L,IC,EXPKJ,TLK
      COMPLEX*16 YLM(0:2*NL_M,-2*NL_M:2*NL_M)
      REAL*8 PI,ATTKJ,XKJ,YKJ,ZKJ,RKJ,ZDKJ,KRKJ,GNT(0:N_GAUNT)
      INTEGER JE,N,JTYP,JNUM,JATL,KTYP,KNUM,KATL,LMJ,LMK
      INTEGER LJ,LK,MJ,MK,L_MIN,L_MAX,L,JLIN,KLIN
      INTEGER IBESS,M
      DATA PI /3.1415926535898D0/
C
      IC=DCMPLX(0.D0,1.D0)
      DO KLIN=1,N
        Y(KLIN)=DCMPLX(0.D0,0.D0)
      ENDDO
      IBESS=3
      JLIN=0
      DO JTYP=1,N_PROT
        LMJ=LMAX(JTYP,JE)
        DO JNUM=1,NATYP(JTYP)
          JATL=NCORR(JNUM,JTYP)
          XJ=SYM_AT(1,JATL)
          YJ=SYM_AT(2,JATL)
          ZJ=SYM_AT(3,JATL)
          DO LJ=0,LMJ
            DO MJ=-LJ,LJ
              JLIN=JLIN+1
              KLIN=0
              DO KTYP=1,N_PROT
                LMK=LMAX(KTYP,JE)
                DO KNUM=1,NATYP(KTYP)
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
                  DO LK=0,LMK
                    L_MIN=ABS(LK-LJ)
                    L_MAX=LK+LJ
                    TLK=DCMPLX(TL(LK,1,KTYP,JE))
                    DO MK=-LK,LK
                      KLIN=KLIN+1
                      IF(KLIN.EQ.JLIN) Y(KLIN)=Y(KLIN)+X(JLIN)
                      IF(KATL.NE.JATL) THEN
                        CALL GAUNT2(LK,MK,LJ,MJ,GNT)
                        SUM_L=DCMPLX(0.D0,0.D0)
                        DO L=L_MIN,L_MAX,2
                          M=MJ-MK
                          IF(ABS(M).LE.L) THEN
                            SUM_L=SUM_L+(IC**L)*HL1(L)*
     1                             YLM(L,M)*GNT(L)
                          ENDIF
                        ENDDO
                        SUM_L=SUM_L*ATTKJ*4.D0*PI*IC
                        Y(KLIN)=Y(KLIN)-TLK*SUM_L*X(JLIN)
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF(JLIN.NE.N) THEN
        WRITE(6,*) 'MS_MATRIX_ACTION dimension mismatch:',JLIN,N
      ENDIF
      RETURN
      END

C=======================================================================
C
      SUBROUTINE MS_MATRIX_ACTION_BLOCK(JE,N,NRHS,X,Y)
C
C   Apply one streamed interaction assembly to several source vectors.
C   Geometry is evaluated once per ordered atom pair. The real Gaunt table
C   is reused while the maximum angular momentum remains unchanged.
C
      USE DIM_MOD
      USE COOR_MOD
      USE INIT_L_MOD
      USE TRANS_MOD
      INTEGER JE,N,NRHS,JTYP,JNUM,JATL,KTYP,KNUM,KATL,LMJ,LMK
      INTEGER LJ,LK,MJ,MK,L_MIN,L_MAX,L,JLIN,KLIN,IBESS,M,IRHS
      INTEGER PROFILE,CLOCK0,CLOCK1,CLOCK_RATE,IPROFLEN,IPROFSTAT
      INTEGER NATOMS,QJ,QK,QSUM,QSELF,NGEOM,NGAUNTCALL
      INTEGER LMAXMAX,QMAX,NGCACHE,JLOC,KLOC,JBASE,KBASE
      INTEGER NGEOMNEW,NGAUNTNEW,ICH
      INTEGER SHELL_OPT,SHELLOPT_STATUS,SHELLOPT_LEN
      INTEGER NTHREAD,THREADREQ,THREADLEN,THREADSTAT,IOS
      INTEGER JATOM,KATOM
      INTEGER OMP_GET_MAX_THREADS
      INTEGER PAIRQ,PAIRL,PAIRLM,PAIRLOC
      INTEGER*8 CACHEBYTES,CACHELIMIT
      COMPLEX*16 X(LINMAX*NATCLU_M,NRHS)
      COMPLEX*16 Y(LINMAX*NATCLU_M,NRHS)
      COMPLEX*16 HL1(0:2*NL_M),SUM_L,IC,EXPKJ,TLK
      COMPLEX*16 YLM(0:2*NL_M,-2*NL_M:2*NL_M)
      REAL*8 PI,ATTKJ,XJ,YJ,ZJ,XKJ,YKJ,ZKJ,RKJ,ZDKJ,KRKJ
      REAL*8 ACTION_SEC
      LOGICAL CACHE_OK
      LOGICAL CACHE_ENABLED
      LOGICAL CACHE_SINGLE
      CHARACTER*1 CACHEVAL
      CHARACTER*1 CACHEPRECVAL
      CHARACTER*1 SHELLOPT_ENV
      INTEGER CACHELEN,CACHESTAT
      LOGICAL AXIAL_ENABLED
      INTEGER AXCHECK,AXVARC,CHKLEN0,CHKSTAT0,AXVLEN0,AXVSTAT0
      CHARACTER*1 CHKVAL0
      CHARACTER*8 AXVSTR0
      COMMON /MSAXENV/ AXCHECK,AXVARC
      CHARACTER*1 AXIALVAL
      INTEGER AXIALLEN,AXIALSTAT
      INTEGER CACHEPRECLEN,CACHEPRECSTAT
      INTEGER, SAVE :: OLDJE=-1,OLDNATOMS=-1,OLDLMAX=-1
      INTEGER, ALLOCATABLE, SAVE :: OLDTYPE(:),OLDNUM(:),OLDLM(:)
      REAL*8, ALLOCATABLE, SAVE :: OLDCOORD(:,:)
      COMPLEX*16, SAVE :: OLDVK=DCMPLX(0.D0,0.D0)
      LOGICAL, SAVE :: GEOMVALID=.FALSE.
      COMPLEX*16, ALLOCATABLE, SAVE :: PAIRY(:,:,:),PAIRH(:,:,:)
      COMPLEX, ALLOCATABLE, SAVE :: PAIRY_SP(:,:,:),PAIRH_SP(:,:,:)
      REAL*8, ALLOCATABLE, SAVE :: PAIRA(:,:)
      CHARACTER*1 PROFVAL
      CHARACTER*32 THREADVAL
      INTEGER, ALLOCATABLE :: LCH(:),MCH(:)
      INTEGER, ALLOCATABLE :: ATOM_TYPE(:),ATOM_NUM(:),ATOM_BASE(:)
      INTEGER, ALLOCATABLE :: ATOM_Q(:),ATOM_LMAX(:)
      REAL*8, ALLOCATABLE, SAVE :: GNTAB(:,:,:)
      REAL*8, ALLOCATABLE :: GNTWORK(:)
      EXTERNAL OMP_GET_MAX_THREADS
      INTEGER, SAVE :: GNTAB_LMAX=-1
      DATA PI /3.1415926535898D0/
C
      PROFILE=0
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_ACTION_PROFILE',
     1     PROFVAL,LENGTH=IPROFLEN,STATUS=IPROFSTAT)
      IF(IPROFSTAT.EQ.0.AND.IPROFLEN.GT.0) THEN
        IF(PROFVAL(1:1).EQ.'1') PROFILE=1
      ENDIF
      THREADREQ=1
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_ACTION_THREADS',
     1     THREADVAL,LENGTH=THREADLEN,STATUS=THREADSTAT)
      IF(THREADSTAT.EQ.0.AND.THREADLEN.GT.0) THEN
        READ(THREADVAL,*,IOSTAT=IOS) THREADREQ
        IF(IOS.NE.0) THREADREQ=1
      ENDIF
      THREADREQ=MAX(1,THREADREQ)
      NTHREAD=MIN(THREADREQ,OMP_GET_MAX_THREADS())
      NTHREAD=MAX(1,NTHREAD)
      NATOMS=0
      LMAXMAX=0
      DO JTYP=1,N_PROT
        LMJ=LMAX(JTYP,JE)
        NATOMS=NATOMS+NATYP(JTYP)
        LMAXMAX=MAX(LMAXMAX,LMJ)
      ENDDO
      QMAX=(LMAXMAX+1)*(LMAXMAX+1)
      NGCACHE=2*LMAXMAX
      IF(PROFILE.EQ.1) THEN
        QSUM=0
        QSELF=0
        DO JTYP=1,N_PROT
          LMJ=LMAX(JTYP,JE)
          QJ=(LMJ+1)*(LMJ+1)
          QSUM=QSUM+NATYP(JTYP)*QJ
          QSELF=QSELF+NATYP(JTYP)*QJ*QJ
        ENDDO
        NGEOM=QSUM*(NATOMS-1)
        NGAUNTCALL=QSUM*QSUM-QSELF
        NGEOMNEW=NATOMS*(NATOMS-1)
        NGAUNTNEW=QMAX*QMAX
        WRITE(6,*) 'TREE_ACTION_PROFILE_BASE_CALLS',NGEOM,
     1             NGAUNTCALL
        WRITE(6,*) 'TREE_ACTION_PROFILE_CACHED_CALLS',NGEOMNEW,
     1             NGAUNTNEW
        WRITE(6,*) 'TREE_ACTION_PROFILE_CACHE_BYTES',
     1             8*(NGCACHE+1)*QMAX*QMAX
      ENDIF
      IF(PROFILE.EQ.1) CALL SYSTEM_CLOCK(CLOCK0,CLOCK_RATE)
C
C     Build angle-independent Gaunt coefficients on a cache miss.
C     The maximum required L is LJ+LK, no larger than 2*LMAXMAX.
      IF(.NOT.ALLOCATED(GNTAB).OR.GNTAB_LMAX.NE.LMAXMAX) THEN
        IF(ALLOCATED(GNTAB)) DEALLOCATE(GNTAB)
        ALLOCATE(LCH(QMAX),MCH(QMAX))
        ICH=0
        DO LJ=0,LMAXMAX
          DO MJ=-LJ,LJ
            ICH=ICH+1
            LCH(ICH)=LJ
            MCH(ICH)=MJ
          ENDDO
        ENDDO
        ALLOCATE(GNTAB(0:NGCACHE,QMAX,QMAX))
        ALLOCATE(GNTWORK(0:N_GAUNT))
        DO JLOC=1,QMAX
          DO KLOC=1,QMAX
            GNTWORK=0.D0
            CALL GAUNT2(LCH(KLOC),MCH(KLOC),LCH(JLOC),MCH(JLOC),
     1                  GNTWORK)
            GNTAB(:,JLOC,KLOC)=GNTWORK(0:NGCACHE)
          ENDDO
        ENDDO
        DEALLOCATE(GNTWORK,LCH,MCH)
        GNTAB_LMAX=LMAXMAX
      ENDIF
      IC=DCMPLX(0.D0,1.D0)
      IBESS=3
      JBASE=0
      JATOM=0
      ALLOCATE(ATOM_TYPE(NATOMS),ATOM_NUM(NATOMS),
     1         ATOM_BASE(NATOMS),ATOM_Q(NATOMS),ATOM_LMAX(NATOMS))
      DO JTYP=1,N_PROT
        LMJ=LMAX(JTYP,JE)
        QJ=(LMJ+1)*(LMJ+1)
        DO JNUM=1,NATYP(JTYP)
          JATOM=JATOM+1
          ATOM_TYPE(JATOM)=JTYP
          ATOM_NUM(JATOM)=JNUM
          ATOM_BASE(JATOM)=JBASE
          ATOM_Q(JATOM)=QJ
          ATOM_LMAX(JATOM)=LMJ
          JBASE=JBASE+QJ
        ENDDO
      ENDDO
      IF(JBASE.NE.N) THEN
        WRITE(6,*) 'MS_MATRIX_ACTION_BLOCK dimension mismatch:',
     1             JBASE,N
      ENDIF
C
C     Reuse exact pair geometry across repeated GMRES operator actions.
C     This cache contains spherical harmonics, spherical Bessel factors,
C     and attenuation; energy-dependent T-matrix factors remain live.
      PAIRLM=2*LMAXMAX
      PAIRQ=(PAIRLM+1)*(PAIRLM+1)
      PAIRL=PAIRLM+1
      CACHELIMIT=536870912
      CACHEBYTES=16_8*INT(NATOMS,8)*INT(NATOMS,8)*
     1           INT(PAIRQ+PAIRL,8)+8_8*INT(NATOMS,8)*
     2           INT(NATOMS,8)
      CACHEBYTES=CACHEBYTES+24_8*INT(NATOMS,8)
      CACHEBYTES=CACHEBYTES+24_8*INT(NATOMS,8)
      CACHE_SINGLE=.FALSE.
      CACHEPRECVAL='0'
      CALL GET_ENVIRONMENT_VARIABLE(
     1     'MSSPEC_TREE_CACHE_SINGLE_PRECISION',CACHEPRECVAL,
     2     LENGTH=CACHEPRECLEN,STATUS=CACHEPRECSTAT)
      IF(CACHEPRECSTAT.EQ.0.AND.CACHEPRECLEN.GT.0) THEN
        IF(CACHEPRECVAL(1:1).EQ.'1') CACHE_SINGLE=.TRUE.
      ENDIF
      IF(CACHE_SINGLE) CACHEBYTES=8_8*INT(NATOMS,8)*
     1     INT(NATOMS,8)*INT(PAIRQ+PAIRL,8)+
     2     8_8*INT(NATOMS,8)*INT(NATOMS,8)+
     3     24_8*INT(NATOMS,8)+24_8*INT(NATOMS,8)
      CACHE_ENABLED=.FALSE.
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_GEOMETRY_CACHE',
     1     CACHEVAL,LENGTH=CACHELEN,STATUS=CACHESTAT)
      IF(CACHESTAT.EQ.0.AND.CACHELEN.GT.0) THEN
        IF(CACHEVAL(1:1).EQ.'1') CACHE_ENABLED=.TRUE.
      ENDIF
C     Exact rotation/axial translation is the default action; set
C     MSSPEC_TREE_AXIAL_TRANSLATION=0 to force the direct Gaunt action.
      AXIAL_ENABLED=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_AXIAL_TRANSLATION',
     1     AXIALVAL,LENGTH=AXIALLEN,STATUS=AXIALSTAT)
      IF(AXIALSTAT.EQ.0.AND.AXIALLEN.GT.0) THEN
        IF(AXIALVAL(1:1).EQ.'0') AXIAL_ENABLED=.FALSE.
      ENDIF
C     Read the pair-action switches once here (serial region), not per
C     pair inside the threaded loop.
      AXCHECK=0
      AXVARC=1
      CHKVAL0='0'
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_AXIAL_VALIDATE',
     1     CHKVAL0,LENGTH=CHKLEN0,STATUS=CHKSTAT0)
      IF(CHKSTAT0.EQ.0.AND.CHKLEN0.GT.0) THEN
        IF(CHKVAL0(1:1).EQ.'1') AXCHECK=1
      ENDIF
      AXVSTR0=' '
      CALL GET_ENVIRONMENT_VARIABLE('MSSPEC_TREE_AXIAL_VARIANT',
     1     AXVSTR0,LENGTH=AXVLEN0,STATUS=AXVSTAT0)
      IF(AXVSTAT0.EQ.0.AND.AXVLEN0.GT.0) THEN
        READ(AXVSTR0,*,IOSTAT=AXVSTAT0) AXVARC
        IF(AXVSTAT0.NE.0) AXVARC=1
      ENDIF
      AXVARC=MAX(0,MIN(3,AXVARC))
C     Diagnostic-only same-build control: use the full-shell loop bounds.
      SHELL_OPT=1
      SHELLOPT_ENV='0'
      CALL GET_ENVIRONMENT_VARIABLE(
     1     'MSSPEC_TREE_AXIAL_FULL_LOOP_REFERENCE',SHELLOPT_ENV,
     2     LENGTH=SHELLOPT_LEN,STATUS=SHELLOPT_STATUS)
      IF(SHELLOPT_STATUS.EQ.0.AND.SHELLOPT_LEN.GT.0) THEN
        IF(SHELLOPT_ENV(1:1).EQ.'1') SHELL_OPT=0
      ENDIF
      CACHE_OK=.FALSE.
      IF(CACHE_ENABLED.AND.CACHEBYTES.LE.CACHELIMIT) THEN
        CACHE_OK=GEOMVALID.AND.OLDJE.EQ.JE.AND.
     1           OLDNATOMS.EQ.NATOMS.AND.OLDLMAX.EQ.LMAXMAX
        IF(CACHE_OK) THEN
          IF(OLDVK.NE.DCMPLX(VK(JE))) CACHE_OK=.FALSE.
          DO JATOM=1,NATOMS
            JTYP=ATOM_TYPE(JATOM)
            JNUM=ATOM_NUM(JATOM)
            JATL=NCORR(JNUM,JTYP)
            IF(OLDTYPE(JATOM).NE.JTYP.OR.OLDNUM(JATOM).NE.JNUM.OR.
     1         OLDLM(JATOM).NE.ATOM_LMAX(JATOM)) CACHE_OK=.FALSE.
            IF(OLDCOORD(1,JATOM).NE.DBLE(SYM_AT(1,JATL)).OR.
     1         OLDCOORD(2,JATOM).NE.DBLE(SYM_AT(2,JATL)).OR.
     2         OLDCOORD(3,JATOM).NE.DBLE(SYM_AT(3,JATL)))
     3         CACHE_OK=.FALSE.
          ENDDO
        ENDIF
        IF(.NOT.CACHE_OK) THEN
          GEOMVALID=.FALSE.
          IF(ALLOCATED(PAIRY)) DEALLOCATE(PAIRY)
          IF(ALLOCATED(PAIRH)) DEALLOCATE(PAIRH)
          IF(ALLOCATED(PAIRY_SP)) DEALLOCATE(PAIRY_SP)
          IF(ALLOCATED(PAIRH_SP)) DEALLOCATE(PAIRH_SP)
          IF(ALLOCATED(PAIRA)) DEALLOCATE(PAIRA)
          IF(ALLOCATED(OLDTYPE)) DEALLOCATE(OLDTYPE)
          IF(ALLOCATED(OLDNUM)) DEALLOCATE(OLDNUM)
          IF(ALLOCATED(OLDLM)) DEALLOCATE(OLDLM)
          IF(ALLOCATED(OLDCOORD)) DEALLOCATE(OLDCOORD)
          IF(CACHE_SINGLE) THEN
            ALLOCATE(PAIRY_SP(0:PAIRQ-1,NATOMS,NATOMS),STAT=IOS)
            IF(IOS.EQ.0) ALLOCATE(
     1           PAIRH_SP(0:PAIRLM,NATOMS,NATOMS),STAT=IOS)
          ELSE
            ALLOCATE(PAIRY(0:PAIRQ-1,NATOMS,NATOMS),STAT=IOS)
            IF(IOS.EQ.0) ALLOCATE(
     1           PAIRH(0:PAIRLM,NATOMS,NATOMS),STAT=IOS)
          ENDIF
          IF(IOS.EQ.0) ALLOCATE(PAIRA(NATOMS,NATOMS),STAT=IOS)
          IF(IOS.EQ.0) ALLOCATE(OLDTYPE(NATOMS),OLDNUM(NATOMS),
     1                          OLDLM(NATOMS),OLDCOORD(3,NATOMS),
     2                          STAT=IOS)
          IF(IOS.EQ.0) THEN
            IF(CACHE_SINGLE) THEN
              PAIRY_SP=CMPLX(0.,0.)
              PAIRH_SP=CMPLX(0.,0.)
            ELSE
              PAIRY=DCMPLX(0.D0,0.D0)
              PAIRH=DCMPLX(0.D0,0.D0)
            ENDIF
            PAIRA=0.D0
            DO JATOM=1,NATOMS
              JTYP=ATOM_TYPE(JATOM)
              JNUM=ATOM_NUM(JATOM)
              JATL=NCORR(JNUM,JTYP)
              OLDTYPE(JATOM)=JTYP
              OLDNUM(JATOM)=JNUM
              OLDLM(JATOM)=ATOM_LMAX(JATOM)
              OLDCOORD(1,JATOM)=DBLE(SYM_AT(1,JATL))
              OLDCOORD(2,JATOM)=DBLE(SYM_AT(2,JATL))
              OLDCOORD(3,JATOM)=DBLE(SYM_AT(3,JATL))
            ENDDO
            DO JATOM=1,NATOMS
              JTYP=ATOM_TYPE(JATOM)
              JNUM=ATOM_NUM(JATOM)
              JATL=NCORR(JNUM,JTYP)
              XJ=SYM_AT(1,JATL)
              YJ=SYM_AT(2,JATL)
              ZJ=SYM_AT(3,JATL)
              LMJ=ATOM_LMAX(JATOM)
              DO KATOM=1,NATOMS
                IF(KATOM.NE.JATOM) THEN
                  KTYP=ATOM_TYPE(KATOM)
                  KNUM=ATOM_NUM(KATOM)
                  KATL=NCORR(KNUM,KTYP)
                  LMK=ATOM_LMAX(KATOM)
                  XKJ=DBLE(SYM_AT(1,KATL)-XJ)
                  YKJ=DBLE(SYM_AT(2,KATL)-YJ)
                  ZKJ=DBLE(SYM_AT(3,KATL)-ZJ)
                  RKJ=DSQRT(XKJ*XKJ+YKJ*YKJ+ZKJ*ZKJ)
                  KRKJ=DBLE(VK(JE))*RKJ
                  ATTKJ=DEXP(-DIMAG(DCMPLX(VK(JE)))*RKJ)
                  EXPKJ=(XKJ+IC*YKJ)/RKJ
                  ZDKJ=ZKJ/RKJ
                  CALL SPH_HAR2(2*NL_M,ZDKJ,EXPKJ,YLM,LMJ+LMK)
                  IBESS=3
                  CALL BESPHE2(LMJ+LMK+1,IBESS,KRKJ,HL1)
                  PAIRA(JATOM,KATOM)=ATTKJ
                  DO L=0,LMJ+LMK
                    IF(CACHE_SINGLE) THEN
                      PAIRH_SP(L,JATOM,KATOM)=HL1(L)
                    ELSE
                      PAIRH(L,JATOM,KATOM)=HL1(L)
                    ENDIF
                    DO M=-L,L
                      PAIRLOC=L*L+L+M
                      IF(CACHE_SINGLE) THEN
                        PAIRY_SP(PAIRLOC,JATOM,KATOM)=YLM(L,M)
                      ELSE
                        PAIRY(PAIRLOC,JATOM,KATOM)=YLM(L,M)
                      ENDIF
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
            OLDJE=JE
            OLDNATOMS=NATOMS
            OLDLMAX=LMAXMAX
            OLDVK=DCMPLX(VK(JE))
            GEOMVALID=.TRUE.
            CACHE_OK=.TRUE.
          ELSE
            IF(ALLOCATED(PAIRY)) DEALLOCATE(PAIRY)
            IF(ALLOCATED(PAIRH)) DEALLOCATE(PAIRH)
            IF(ALLOCATED(PAIRY_SP)) DEALLOCATE(PAIRY_SP)
            IF(ALLOCATED(PAIRH_SP)) DEALLOCATE(PAIRH_SP)
            IF(ALLOCATED(PAIRA)) DEALLOCATE(PAIRA)
            IF(ALLOCATED(OLDTYPE)) DEALLOCATE(OLDTYPE)
            IF(ALLOCATED(OLDNUM)) DEALLOCATE(OLDNUM)
            IF(ALLOCATED(OLDLM)) DEALLOCATE(OLDLM)
            IF(ALLOCATED(OLDCOORD)) DEALLOCATE(OLDCOORD)
          ENDIF
        ENDIF
      ELSE
        GEOMVALID=.FALSE.
        IF(ALLOCATED(PAIRY)) DEALLOCATE(PAIRY)
        IF(ALLOCATED(PAIRH)) DEALLOCATE(PAIRH)
        IF(ALLOCATED(PAIRY_SP)) DEALLOCATE(PAIRY_SP)
        IF(ALLOCATED(PAIRH_SP)) DEALLOCATE(PAIRH_SP)
        IF(ALLOCATED(PAIRA)) DEALLOCATE(PAIRA)
        IF(ALLOCATED(OLDTYPE)) DEALLOCATE(OLDTYPE)
        IF(ALLOCATED(OLDNUM)) DEALLOCATE(OLDNUM)
        IF(ALLOCATED(OLDLM)) DEALLOCATE(OLDLM)
        IF(ALLOCATED(OLDCOORD)) DEALLOCATE(OLDCOORD)
      ENDIF
      IF(PROFILE.EQ.1) WRITE(6,*)
     1  'TREE_ACTION_GEOMETRY_CACHE',CACHE_OK,CACHE_ENABLED,
     2  CACHEBYTES
      IF(PROFILE.EQ.1) WRITE(6,*)
     1  'TREE_ACTION_CACHE_SINGLE_PRECISION',CACHE_SINGLE
C
C     Give each thread complete target-atom rows.  This removes the
C     replicated N*NRHS*NTHREAD output slabs and the final reduction.
C     Source atoms are summed in their original order for each target.
      Y(1:N,1:NRHS)=DCMPLX(0.D0,0.D0)
C$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(NTHREAD)
C$OMP& PRIVATE(JATOM,KATOM,JTYP,JNUM,JATL,KTYP,KNUM,KATL,LMJ,LMK,
C$OMP& PAIRLOC,QJ,QK,JBASE,KBASE,JLOC,KLOC,LJ,LK,MJ,MK,L_MIN,L_MAX,
C$OMP& L,M,JLIN,KLIN,IRHS,IBESS,TLK,SUM_L,XJ,YJ,ZJ,XKJ,YKJ,ZKJ,
C$OMP& RKJ,KRKJ,ATTKJ,EXPKJ,ZDKJ,YLM,HL1)
C$OMP DO SCHEDULE(STATIC)
      DO KATOM=1,NATOMS
        KTYP=ATOM_TYPE(KATOM)
        KNUM=ATOM_NUM(KATOM)
        KBASE=ATOM_BASE(KATOM)
        QK=ATOM_Q(KATOM)
        LMK=ATOM_LMAX(KATOM)
        KATL=NCORR(KNUM,KTYP)
        DO JATOM=1,NATOMS
          JTYP=ATOM_TYPE(JATOM)
          JNUM=ATOM_NUM(JATOM)
          JBASE=ATOM_BASE(JATOM)
          QJ=ATOM_Q(JATOM)
          LMJ=ATOM_LMAX(JATOM)
          JATL=NCORR(JNUM,JTYP)
          XJ=SYM_AT(1,JATL)
          YJ=SYM_AT(2,JATL)
          ZJ=SYM_AT(3,JATL)
          IF(KATL.EQ.JATL) THEN
            DO JLOC=1,QJ
              JLIN=JBASE+JLOC
              KLIN=KBASE+JLOC
              DO IRHS=1,NRHS
                Y(KLIN,IRHS)=Y(KLIN,IRHS)+X(JLIN,IRHS)
              ENDDO
            ENDDO
          ELSE
            IF(CACHE_OK) THEN
              ATTKJ=PAIRA(JATOM,KATOM)
              DO L=0,LMJ+LMK
                IF(CACHE_SINGLE) THEN
                  HL1(L)=PAIRH_SP(L,JATOM,KATOM)
                ELSE
                  HL1(L)=PAIRH(L,JATOM,KATOM)
                ENDIF
                DO M=-L,L
                  PAIRLOC=L*L+L+M
                  IF(CACHE_SINGLE) THEN
                    YLM(L,M)=PAIRY_SP(PAIRLOC,JATOM,KATOM)
                  ELSE
                    YLM(L,M)=PAIRY(PAIRLOC,JATOM,KATOM)
                  ENDIF
                ENDDO
              ENDDO
            ELSE
              XKJ=DBLE(SYM_AT(1,KATL)-XJ)
              YKJ=DBLE(SYM_AT(2,KATL)-YJ)
              ZKJ=DBLE(SYM_AT(3,KATL)-ZJ)
              RKJ=DSQRT(XKJ*XKJ+YKJ*YKJ+ZKJ*ZKJ)
              KRKJ=DBLE(VK(JE))*RKJ
              ATTKJ=DEXP(-DIMAG(DCMPLX(VK(JE)))*RKJ)
              EXPKJ=(XKJ+IC*YKJ)/RKJ
              ZDKJ=ZKJ/RKJ
              CALL SPH_HAR2(2*NL_M,ZDKJ,EXPKJ,YLM,LMJ+LMK)
              IBESS=3
              CALL BESPHE2(LMJ+LMK+1,IBESS,KRKJ,HL1)
            ENDIF
            IF(AXIAL_ENABLED) THEN
              IF(CACHE_OK) THEN
                XKJ=DBLE(SYM_AT(1,KATL)-XJ)
                YKJ=DBLE(SYM_AT(2,KATL)-YJ)
                ZKJ=DBLE(SYM_AT(3,KATL)-ZJ)
                RKJ=DSQRT(XKJ*XKJ+YKJ*YKJ+ZKJ*ZKJ)
                KRKJ=DBLE(VK(JE))*RKJ
                ATTKJ=PAIRA(JATOM,KATOM)
              ENDIF
              CALL MS_AXIAL_PAIR_ACTION(JE,NRHS,LMJ,LMK,KTYP,
     1             JBASE,KBASE,QMAX,NGCACHE,GNTAB,X,Y,KRKJ,
     2             ATTKJ,XKJ,YKJ,ZKJ,SHELL_OPT)
            ELSE
            JLOC=0
            DO LJ=0,LMJ
              DO MJ=-LJ,LJ
                JLOC=JLOC+1
                JLIN=JBASE+JLOC
                KLOC=0
                DO LK=0,LMK
                  L_MIN=ABS(LK-LJ)
                  L_MAX=LK+LJ
                  TLK=DCMPLX(TL(LK,1,KTYP,JE))
                  DO MK=-LK,LK
                    KLOC=KLOC+1
                    KLIN=KBASE+KLOC
                    SUM_L=DCMPLX(0.D0,0.D0)
                    DO L=L_MIN,L_MAX,2
                      M=MJ-MK
                      IF(ABS(M).LE.L) THEN
                        SUM_L=SUM_L+(IC**L)*HL1(L)*
     1                       YLM(L,M)*GNTAB(L,JLOC,KLOC)
                      ENDIF
                    ENDDO
                    SUM_L=-TLK*SUM_L*ATTKJ*4.D0*PI*IC
                    DO IRHS=1,NRHS
                      Y(KLIN,IRHS)=Y(KLIN,IRHS)+
     1                    SUM_L*X(JLIN,IRHS)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C$OMP END DO
C$OMP END PARALLEL
      DEALLOCATE(ATOM_TYPE,ATOM_NUM,ATOM_BASE,ATOM_Q,ATOM_LMAX)
      IF(PROFILE.EQ.1) THEN
        CALL SYSTEM_CLOCK(CLOCK1)
        ACTION_SEC=DBLE(CLOCK1-CLOCK0)/DBLE(CLOCK_RATE)
        WRITE(6,*) 'TREE_ACTION_PROFILE_SECONDS',ACTION_SEC
      ENDIF
      RETURN
      END

C=======================================================================
C
      SUBROUTINE MS_AXIAL_PAIR_ACTION(JE,NRHS,LMS,LMT,KTYP,
     1                                JBASE,KBASE,QMAX,NGCACHE,
     2                                GNTAB,X,Y,KR,ATT,DX,DY,DZ,
     3                                SHELL_ONLY)
C
C   Exact rotation--axial-translation--inverse-rotation action for one
C   directed spherical, non-spin pair.  The direct path remains default.
C   This prototype uses the repository's DJMN2 Wigner-d recurrence.
C
      USE DIM_MOD
      USE TRANS_MOD
      INTEGER JE,NRHS,LMS,LMT,KTYP,JBASE,KBASE,QMAX,NGCACHE
      INTEGER SHELL_ONLY
      INTEGER QS,QT,LS,LT,MS,MT,L,M,JS,KT,IS,IT,IRHS
      INTEGER JOFF,KOFF,JIDX,KIDX,LMAXROT,AXVAR,AXSTAT,AXLEN
      COMPLEX*16 X(LINMAX*NATCLU_M,NRHS)
      COMPLEX*16 Y(LINMAX*NATCLU_M,NRHS)
      REAL*8 GNTAB(0:NGCACHE,QMAX,QMAX)
      REAL*8 KR,ATT,DX,DY,DZ,THETA,PHI,PI
      REAL RLM(-NL_M+1:NL_M-1,-NL_M+1:NL_M-1,0:NL_M-1)
      COMPLEX*16, ALLOCATABLE :: DS(:,:),DT(:,:),XR(:,:),YR(:,:)
      COMPLEX*16, ALLOCATABLE :: YTMP(:,:),YDIR(:,:)
      COMPLEX*16 H(0:2*NL_M)
      COMPLEX*16 YLM(0:2*NL_M,-2*NL_M:2*NL_M),EXPKJ
      COMPLEX*16 PHASE,COEF,TLK,IC,ZPHI
      COMPLEX*16 DIFFN,DENOM,SUM_L
      LOGICAL CHECKAXIAL
      CHARACTER*1 CHECKVAL
      CHARACTER*8 AXVARSTR
      INTEGER AXCHECK,AXVARC
      COMMON /MSAXENV/ AXCHECK,AXVARC
      EXTERNAL DJMN2
      DATA PI /3.14159265358979323846D0/
      DATA IC /(0.D0,1.D0)/
C
      QS=(LMS+1)*(LMS+1)
      QT=(LMT+1)*(LMT+1)
      CHECKAXIAL=(AXCHECK.EQ.1)
      IF(CHECKAXIAL.AND.KBASE.EQ.0.AND.JBASE.GT.0.AND.
     1   JBASE.LT.2*QMAX) ALLOCATE(YTMP(QT,NRHS),YDIR(QT,NRHS))
      AXVAR=AXVARC
      LMAXROT=MAX(LMS,LMT)
      THETA=DACOS(MAX(-1.D0,MIN(1.D0,DZ/
     1       DSQRT(DX*DX+DY*DY+DZ*DZ))))
      PHI=DATAN2(DY,DX)
      ALLOCATE(DS(QS,QS),DT(QT,QT),XR(QS,NRHS),YR(QT,NRHS))
      DS=DCMPLX(0.D0,0.D0)
      DT=DCMPLX(0.D0,0.D0)
C
C     D_l(m_out,m_in)=exp(+i*m_out*phi) d_l(m_out,m_in,+theta).
C     This is the convention verified against the direct Gaunt expansion.
      DO LS=0,LMS
        IF(MOD(AXVAR,2).EQ.0) THEN
          CALL DJMN2(REAL(-THETA),RLM,LS,0)
        ELSE
          CALL DJMN2(REAL(THETA),RLM,LS,0)
        ENDIF
        JOFF=LS*LS
        DO MS=-LS,LS
          IF(AXVAR.LT.2) THEN
            ZPHI=EXP(DCMPLX(0.D0,DBLE(MS)*PHI))
          ELSE
            ZPHI=EXP(DCMPLX(0.D0,-DBLE(MS)*PHI))
          ENDIF
          DO MT=-LS,LS
            JS=JOFF+MS+LS+1
            KT=JOFF+MT+LS+1
            DS(JS,KT)=ZPHI*DBLE(RLM(MT,MS,LS))
          ENDDO
        ENDDO
      ENDDO
      DO LT=0,LMT
        IF(MOD(AXVAR,2).EQ.0) THEN
          CALL DJMN2(REAL(-THETA),RLM,LT,0)
        ELSE
          CALL DJMN2(REAL(THETA),RLM,LT,0)
        ENDIF
        KOFF=LT*LT
        DO MS=-LT,LT
          IF(AXVAR.LT.2) THEN
            ZPHI=EXP(DCMPLX(0.D0,DBLE(MS)*PHI))
          ELSE
            ZPHI=EXP(DCMPLX(0.D0,-DBLE(MS)*PHI))
          ENDIF
          DO MT=-LT,LT
            JS=KOFF+MS+LT+1
            KT=KOFF+MT+LT+1
            DT(JS,KT)=ZPHI*DBLE(RLM(MT,MS,LT))
          ENDDO
        ENDDO
      ENDDO
C
C     Rotate source vector with D_source^T.
      IF(SHELL_ONLY.EQ.1) THEN
        DO IRHS=1,NRHS
          DO LS=0,LMS
            JOFF=LS*LS
            DO MS=-LS,LS
              JS=JOFF+MS+LS+1
              XR(JS,IRHS)=DCMPLX(0.D0,0.D0)
              DO MT=-LS,LS
                IS=JOFF+MT+LS+1
                XR(JS,IRHS)=XR(JS,IRHS)+DS(IS,JS)*
     1                        X(JBASE+IS,IRHS)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO IRHS=1,NRHS
          DO JS=1,QS
            XR(JS,IRHS)=DCMPLX(0.D0,0.D0)
            DO IS=1,QS
              XR(JS,IRHS)=XR(JS,IRHS)+DS(IS,JS)*
     1                    X(JBASE+IS,IRHS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
C     Axial pair block: only equal m couples for z-directed translation.
      CALL BESPHE2(LMS+LMT+1,3,KR,H)
      YR=DCMPLX(0.D0,0.D0)
      IF(SHELL_ONLY.EQ.1) THEN
        DO LT=0,LMT
          TLK=DCMPLX(TL(LT,1,KTYP,JE))
          KOFF=LT*LT
          DO MT=-LT,LT
            KT=KOFF+MT+LT+1
            DO LS=0,LMS
              JOFF=LS*LS
              IF(ABS(MT).LE.LS) THEN
                JS=JOFF+MT+LS+1
                COEF=DCMPLX(0.D0,0.D0)
                DO L=ABS(LS-LT),LS+LT,2
                  COEF=COEF+(IC**L)*H(L)*
     1              DSQRT(DBLE(2*L+1)/(4.D0*PI))*
     2              GNTAB(L,JS,KT)
                ENDDO
                COEF=-TLK*COEF*ATT*4.D0*PI*IC
                DO IRHS=1,NRHS
                  YR(KT,IRHS)=YR(KT,IRHS)+COEF*XR(JS,IRHS)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO LT=0,LMT
          TLK=DCMPLX(TL(LT,1,KTYP,JE))
          KOFF=LT*LT
          DO MT=-LT,LT
            KT=KOFF+MT+LT+1
            DO LS=0,LMS
              JOFF=LS*LS
              DO MS=-LS,LS
                IF(MS.EQ.MT) THEN
                  JS=JOFF+MS+LS+1
                  COEF=DCMPLX(0.D0,0.D0)
                  DO L=ABS(LS-LT),LS+LT,2
                    COEF=COEF+(IC**L)*H(L)*
     1                DSQRT(DBLE(2*L+1)/(4.D0*PI))*
     2                GNTAB(L,JS,KT)
                  ENDDO
                  COEF=-TLK*COEF*ATT*4.D0*PI*IC
                  DO IRHS=1,NRHS
                    YR(KT,IRHS)=YR(KT,IRHS)+COEF*XR(JS,IRHS)
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
C     Rotate target vector back with conjugate(D_target).
      IF(ALLOCATED(YTMP)) THEN
        DO IRHS=1,NRHS
          DO LT=0,LMT
            KOFF=LT*LT
            DO MT=-LT,LT
              KT=KOFF+MT+LT+1
              YTMP(KT,IRHS)=DCMPLX(0.D0,0.D0)
              DO MS=-LT,LT
                JS=KOFF+MS+LT+1
                YTMP(KT,IRHS)=YTMP(KT,IRHS)+
     1                      DCONJG(DT(KT,JS))*YR(JS,IRHS)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        YDIR=DCMPLX(0.D0,0.D0)
        EXPKJ=DCMPLX(DX,DY)/DSQRT(DX*DX+DY*DY+DZ*DZ)
        CALL SPH_HAR2(2*NL_M,DZ/DSQRT(DX*DX+DY*DY+DZ*DZ),
     1                EXPKJ,YLM,LMS+LMT)
        DO LT=0,LMT
          TLK=DCMPLX(TL(LT,1,KTYP,JE))
          KOFF=LT*LT
          DO MT=-LT,LT
            KT=KOFF+MT+LT+1
            DO LS=0,LMS
              JOFF=LS*LS
              DO MS=-LS,LS
                JS=JOFF+MS+LS+1
                SUM_L=DCMPLX(0.D0,0.D0)
                DO L=ABS(LS-LT),LS+LT,2
                  M=MS-MT
                  IF(ABS(M).LE.L) SUM_L=SUM_L+(IC**L)*H(L)*
     1                 YLM(L,M)*GNTAB(L,JS,KT)
                ENDDO
                SUM_L=-TLK*SUM_L*ATT*4.D0*PI*IC
                DO IRHS=1,NRHS
                  YDIR(KT,IRHS)=YDIR(KT,IRHS)+
     1                          SUM_L*X(JBASE+JS,IRHS)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DIFFN=DCMPLX(0.D0,0.D0)
        DENOM=DCMPLX(0.D0,0.D0)
        DO IRHS=1,NRHS
          DO KT=1,QT
            DIFFN=DIFFN+DCONJG(YTMP(KT,IRHS)-YDIR(KT,IRHS))*
     1                     (YTMP(KT,IRHS)-YDIR(KT,IRHS))
            DENOM=DENOM+DCONJG(YDIR(KT,IRHS))*YDIR(KT,IRHS)
          ENDDO
        ENDDO
        WRITE(6,*) 'AXIAL_PAIR_RELERR',DSQRT(DBLE(DIFFN)) /
     1       MAX(DSQRT(DBLE(DENOM)),1.D-300),LMS,LMT,KR
        DEALLOCATE(YTMP,YDIR)
      ENDIF
      IF(SHELL_ONLY.EQ.1) THEN
        DO IRHS=1,NRHS
          DO LT=0,LMT
            KOFF=LT*LT
            DO MT=-LT,LT
              KT=KOFF+MT+LT+1
              COEF=DCMPLX(0.D0,0.D0)
              DO MS=-LT,LT
                JS=KOFF+MS+LT+1
                COEF=COEF+DCONJG(DT(KT,JS))*YR(JS,IRHS)
              ENDDO
              Y(KBASE+KT,IRHS)=Y(KBASE+KT,IRHS)+COEF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO IRHS=1,NRHS
          DO KT=1,QT
            COEF=DCMPLX(0.D0,0.D0)
            DO JS=1,QT
              COEF=COEF+DCONJG(DT(KT,JS))*YR(JS,IRHS)
            ENDDO
            Y(KBASE+KT,IRHS)=Y(KBASE+KT,IRHS)+COEF
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(DS,DT,XR,YR)
      RETURN
      END

