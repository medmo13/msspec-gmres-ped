      SUBROUTINE AMAS(NIV,ATOME,COORD,VALZ,ISURF,COUPUR,ROT,IRE,NATYP,NB
     &Z,NAT2,NCOUCH,NMAX)
C
C  This routine generates a cluster from the knowledge of its
C                          lattice vectors
C
      USE DIM_MOD
C
      USE ADSORB_MOD , NCOUCH1 => NCOUCH
      USE BASES_MOD
      USE MILLER_MOD , IM1 => IH, IM2 => IK, IM3 => II, IM4 => IL
      USE OUTUNITS_MOD
      USE RESEAU_MOD
C
      DIMENSION VALZ(NATCLU_M)
      DIMENSION ROT(3,3),IRE(NATCLU_M,2),NATYP(NATM),ITA(NATCLU_M)
      DIMENSION ATOME(3,NATCLU_M),ATRSU(3,NATCLU_M),COORD(3,NATCLU_M)
      DIMENSION ROTINV(3,3),XINIT(3,1),XFIN(3,1)
C
C
C
      NCOUCH=0
      WRITE(IUO1,10) ISURF
  10  FORMAT(//,18X,'ATOM (0,0,0) ON THE SURFACE PLANE IS OF TYPE ',I2)
      NBZ=0
      CALL INVMAT(ROT,ROTINV)
      IF(IVG0.EQ.0) THEN
        CALL CHBASE(NATP_M,ATBAS)
      ENDIF
      NB1=0
      NB2=0
      DO NTYP=1,NAT
        NBAT=0
        DO NUM=1,NMAX
          NB1=NB1+1
          IRE(NB1,1)=0
          IRE(NB1,2)=0
          IF(IVG0.LE.1) THEN
            CALL NUMAT(NUM,NIV,IA,IB,IC)
          ELSE
            BSURA=1.
            CSURA=1.
          ENDIF
          IF(IVG0.LE.1) THEN
            XA=FLOAT(IA)
            XB=FLOAT(IB)
            XC=FLOAT(IC)
          ELSEIF(IVG0.EQ.2) THEN
            XA=FLOAT(NUM-1)
            XB=FLOAT(NUM-1)
            XC=FLOAT(NUM-1)
          ENDIF
          IF(IVG0.EQ.1) THEN
            IF(IVN(1).EQ.0) THEN
              ITA(NUM)=IA
            ELSEIF(IVN(2).EQ.0) THEN
              ITA(NUM)=IB
            ELSEIF(IVN(3).EQ.0) THEN
              ITA(NUM)=IC
            ENDIF
            IF((ITA(NUM).EQ.ITA(NUM-1)).AND.(NUM.GT.1)) GOTO 30
          ENDIF
          DO J=1,3
            K=J+3*(NTYP-1)
            O=ATBAS(K)
            ATOME(J,NB1)=O+XA*VECBAS(J)+XB*VECBAS(J+3)+XC*VECBAS(J+6)
          ENDDO
          DO I=1,3
            M=I+3*(ISURF-1)
            XINIT(I,1)=ATOME(I,NB1)-ATBAS(M)
          ENDDO
          CALL MULMAT(ROTINV,3,3,XINIT,3,1,XFIN)
          DO I=1,3
            ATRSU(I,NB1)=XFIN(I,1)
          ENDDO
          CALL TEST1(COUPUR,NB1,NB2,ATRSU,COORD,VALZ,NBAT,IRE,NBZ)
  30      CONTINUE
          ENDDO
        NATYP(NTYP)=NBAT
      ENDDO
      IF(IADS.GE.1) THEN
        N0=NBZ
        DO JADS=1,NADS1
          NB1=NB1+1
          DO I=1,3
            COORD(I,NB1)=ADS(I,JADS)
          ENDDO
          N1=0
          DO N=1,N0
            D=ABS(COORD(3,NB1)-VALZ(N))
            IF(D.LT.0.0001) N1=N1+1
          ENDDO
          IF(N1.EQ.0) THEN
            N0=N0+1
            VALZ(N0)=COORD(3,NB1)
          ENDIF
        ENDDO
        NANEW1=NADS1+NADS2
        NATYP(NAT+1)=NADS1
        IF(NANEW1.EQ.NADS1) GOTO 99
        DO JADS=NADS1+1,NANEW1
          NB1=NB1+1
          DO I=1,3
            COORD(I,NB1)=ADS(I,JADS)
          ENDDO
          N1=0
          DO N=1,N0
            D=ABS(COORD(3,NB1)-VALZ(N))
            IF(D.LT.0.0001) N1=N1+1
          ENDDO
          IF(N1.EQ.0) THEN
            N0=N0+1
            VALZ(N0)=COORD(3,NB1)
          ENDIF
        ENDDO
        NATYP(NAT+2)=NADS2
        NANEW2=NANEW1+NADS3
        IF(NANEW2.EQ.NANEW1) GOTO 99
        DO JADS=NANEW1+1,NANEW2
          NB1=NB1+1
          DO I=1,3
            COORD(I,NB1)=ADS(I,JADS)
          ENDDO
          N1=0
          DO N=1,N0
            D=ABS(COORD(3,NB1)-VALZ(N))
            IF(D.LT.0.0001) N1=N1+1
          ENDDO
          IF(N1.EQ.0) THEN
            N0=N0+1
            VALZ(N0)=COORD(3,NB1)
          ENDIF
        ENDDO
        NATYP(NAT+3)=NADS3
  99    CONTINUE
        NCOUCH=N0-NBZ
        NBZ=N0
      ENDIF
C
      RETURN
C
      END
