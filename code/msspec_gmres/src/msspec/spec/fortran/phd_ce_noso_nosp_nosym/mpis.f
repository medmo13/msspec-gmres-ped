C
C
C======================================================================
C
      SUBROUTINE MPIS(N,NLM,ITYP,IGS,JE,QI,TAU)
C
C
C   This subroutine construct the correlation matrices and uses
C    LU decomposition method to do the matrix inversion.
C   The inverse matrix which is the contribution of a small atom group 
C    is kept for further use.
C
C                                                  H. -F. Zhao : 2007
C
C                                    Last modified (DS) : 13 May 2009
C
      USE DIM_MOD
      USE COOR_MOD
      USE INIT_L_MOD
      USE GAUNT_C_MOD
      USE TRANS_MOD
      USE CORREXP_MOD
C
      INTEGER NLM(NGR_M),ITYP(NGR_M),IGS(NGR_M)
      COMPLEX*16 TAU(LINMAX,LINFMAX,NATCLU_M)
C
      REAL QI
C
      COMPLEX*16 ZEROC,ONEC,IC
C
      COMPLEX*16 ATTL(0:NT_M,NATM)
      COMPLEX*16 EXPJN,ATTJN
      COMPLEX*16 YLM(0:NLTWO,-NLTWO:NLTWO)
      COMPLEX*16 HL1(0:NLTWO)
      COMPLEX*16 SUM_L,SUM_L2
      COMPLEX*16 SUM_L_A,SUM_L2_A,SUM_L_B,SUM_L2_B
C
      REAL*8 FOURPI
      REAL*8 XJN,YJN,ZJN,RJN,KRJN,ZDJN
      REAL*8 IM_VK,RE_VK
C
      INTEGER IPIV(NLMM),ONE_L,IN1
C
      COMPLEX*16 FOURPI_IC,IC_L,IC_REF,TEMP,TEMP1,TEMP2,CN1
      COMPLEX*16 AINV(NLMM,NLMM),IN(NLMM,LINFMAX)
C
      DATA FOURPI /12.566370614359D0/
C
      ZEROC=(0.D0,0.D0)
      ONEC=(1.D0,0.D0)
      IC=(0.D0,1.D0)
      IBESS=3
      FOURPI_IC=-IC*FOURPI
C
      LM0=LMAX(1,JE)
      LM0=MIN(LM0,LF2)
      NRHS=(LM0+1)*(LM0+1)
      INDJ=0
C
      NM=0
      DO I=1,N-1
        J=NLM(I)+1
        NM=NM+J*J
      ENDDO
      L=NLM(N)
      LNMAX=L
      L=(L+1)*(L+1)
      NM1=NM+1
      NML=NM+L
      NTYP=ITYP(N)
C
      DO L=0,LNMAX
        ATTL(L,N)=DCMPLX(TL(L,1,NTYP,JE))
      ENDDO
      IM_VK=-DIMAG(DCMPLX(VK(JE)))
      RE_VK=DBLE(VK(JE))
C
C     set up matrix blocks C((N-1)*1) and D(1*(N-1))
C
      I=IGS(N)
      XN=SYM_AT(1,I)
      YN=SYM_AT(2,I)
      ZN=SYM_AT(3,I)
C
      DO J=1,N-1
        JATL=IGS(J)
        LJMAX=NLM(J)
        JTYP=ITYP(J)
        J1=J-1
C
        XJN=DBLE(SYM_AT(1,JATL)-XN)
        YJN=DBLE(SYM_AT(2,JATL)-YN)
        ZJN=DBLE(SYM_AT(3,JATL)-ZN)
        RJN=DSQRT(XJN*XJN+YJN*YJN+ZJN*ZJN)
        KRJN=RE_VK*RJN
        ATTJN=FOURPI_IC*DEXP(IM_VK*RJN)
        EXPJN=(XJN+IC*YJN)/RJN
        ZDJN=ZJN/RJN
        CALL SPH_HAR2(2*NL_M,ZDJN,EXPJN,YLM,LNMAX+LJMAX)
        CALL BESPHE2(LNMAX+LJMAX+1,IBESS,KRJN,HL1)
        DO L=0,LJMAX
          ATTL(L,J)=ATTJN*DCMPLX(TL(L,1,JTYP,JE))
        ENDDO
C
        II=NM
        IN1=-1
        CN1=IC
        JJ=0
C
        DO LN=0,LNMAX
          ILN=LN*LN+LN+1
          IN1=-IN1
          CN1=-CN1*IC
C
          DO MLN=-LN,LN
            INDN=ILN+MLN
            II=II+1
            JJ0=J1*INDJ
            ONE_L=-IN1
            IC_REF=-CN1*IC
C
            DO LJ=0,LJMAX
              ILJ=LJ*LJ+LJ+1
              L_MIN=ABS(LJ-LN)
              L_MAX=LJ+LN
              ONE_L=-ONE_L
              IC_REF=IC_REF*IC
C
C  Case MLJ equal to zero
C 
              JJ1=JJ0+ILJ
              IF(LJ.GE.LN) THEN
                IC_L=-IC_REF
              ELSE
                IC_L=-ONEC/IC_REF
              ENDIF
C
              SUM_L=ZEROC
              SUM_L2=ZEROC
C
              DO L=L_MIN,L_MAX,2
                IC_L=-IC_L
                IF(ABS(MLN).LE.L) THEN
                  TEMP=IC_L*HL1(L)*GNT(L,ILJ,INDN)
                  SUM_L=SUM_L+YLM(L,MLN)*TEMP
                  SUM_L2=SUM_L2+DCONJG(YLM(L,MLN))*TEMP
                ENDIF
              ENDDO
C
              IF(ONE_L.EQ.-1) SUM_L2=-SUM_L2
              A(JJ1,II)=ATTL(LJ,J)*SUM_L
              A(II,JJ1)=ATTJN*ATTL(LN,N)*SUM_L2
C
C
C  Case MLJ not equal to zero
C 
              DO MLJ=1,LJ
                INDJ=ILJ+MLJ
                INDJN=ILJ-MLJ
                JJ1=JJ0+INDJ
                JJ1N=JJ0+INDJN
                MA=MLN-MLJ
                MB=MLN+MLJ
                IF(LJ.GE.LN) THEN
                  IC_L=-IC_REF
                ELSE
                  IC_L=-ONEC/IC_REF
                ENDIF
C
                SUM_L_A=ZEROC
                SUM_L2_A=ZEROC
                SUM_L_B=ZEROC
                SUM_L2_B=ZEROC
C
                DO L=L_MIN,L_MAX,2
                  IC_L=-IC_L
                  IF(ABS(MA).LE.L) THEN
                    TEMP1=IC_L*HL1(L)*GNT(L,INDJ,INDN)
                    SUM_L_A=SUM_L_A+YLM(L,MA)*TEMP1
                    SUM_L2_A=SUM_L2_A+DCONJG(YLM(L,MA))*TEMP1
                  ENDIF
                  IF(ABS(MB).LE.L) THEN
                    TEMP2=IC_L*HL1(L)*GNT(L,INDJN,INDN)
                    SUM_L_B=SUM_L_B+YLM(L,MB)*TEMP2
                    SUM_L2_B=SUM_L2_B+DCONJG(YLM(L,MB))*TEMP2
                  ENDIF
                ENDDO
C
                IF(ONE_L.EQ.-1) THEN
                   SUM_L2_A=-SUM_L2_A
                   SUM_L2_B=-SUM_L2_B
                ENDIF
                A(JJ1,II)=ATTL(LJ,J)*SUM_L_A
                A(II,JJ1)=ATTJN*ATTL(LN,N)*SUM_L2_A
                A(JJ1N,II)=ATTL(LJ,J)*SUM_L_B
                A(II,JJ1N)=ATTJN*ATTL(LN,N)*SUM_L2_B
              ENDDO
C
C
            ENDDO
            JJ=JJ0+INDJ
C
          ENDDO
        ENDDO
C
        JJ=JJ-INDN
C
      ENDDO
C
C     add B to A
C
      DO I=NM1,NML
        DO J=NM1,NML
          IF(J.EQ.I) THEN
            A(J,I)=ONEC
          ELSE
            A(J,I)=ZEROC
          ENDIF
        ENDDO
      ENDDO
C
C     construct AINV
C
      DO I=1,NML
        DO J=1,NML
          AINV(J,I)=A(J,I)
        ENDDO
      ENDDO
C
C
C          matrix inversion(ax=b)
C
      CALL ZGETRF(NML,NML,AINV,NLMM,IPIV,INFO1)
      IF(INFO1.NE.0) THEN
        WRITE(6,*) '     --->  INFO1 =',INFO1
      ELSE
C
        DO I=1,NRHS
          DO J=1,NML
            IF(J.EQ.I) THEN
              IN(J,I)=(1.D0,0.D0)
            ELSE
              IN(J,I)=(0.D0,0.D0)
            ENDIF
          ENDDO
        ENDDO
C
        CALL ZGETRS('N',NML,NRHS,AINV,NLMM,IPIV,IN,NLMM,INFO)
        IF(INFO.NE.0) THEN
          WRITE(6,*) '    --->  INFO =',INFO
        ENDIF
      ENDIF
C
C      sum of tau
C
      KLIN=0
      DO K=1,N
        KATL=IGS(K)
        LMK=NLM(K)
        INDKM=(LMK+1)*(LMK+1)
C
        DO INDJ=1,NRHS
C
          DO INDK=1,INDKM
            KLIN=KLIN+1
C
            TAU(INDK,INDJ,KATL)=TAU(INDK,INDJ,KATL)
     1                 +DBLE(QI)*IN(KLIN,INDJ)
C
          ENDDO
          KLIN=KLIN-INDKM
C
        ENDDO
        KLIN=KLIN+INDKM
C
      ENDDO
C
      RETURN
C
      END
