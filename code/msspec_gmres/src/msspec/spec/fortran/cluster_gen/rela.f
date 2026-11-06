C
C=======================================================================
C
      SUBROUTINE RELA(NINI,NFIN,NAT,VALINI,VALFIN,VALFIN2,COORD,NTYP,REL
     &,L)
C
      USE DIM_MOD
C
      USE ADSORB_MOD , I1 => IADS, N1 => NADS1, N2 => NADS2, N3 => NADS3
     &
      USE OUTUNITS_MOD
      USE RELADS_MOD
      USE RELAX_MOD
C
      DIMENSION VALINI(NATCLU_M),VALFIN(NATCLU_M),REL(NATCLU_M)
      DIMENSION NTYP(NATM),COORD(3,NATCLU_M),LSP(2),DZA(2),DZB(2)
      DIMENSION DYA(2),DYB(2),VALFIN2(NATCLU_M),KZ(1000)
C
      DATA SMALL /0.0001/
C
      IF((IREL.EQ.1).OR.((IREL.EQ.0).AND.(NRELA.GT.0))) THEN
C
        CALL ORDRE(NINI,VALINI,NFIN,VALFIN)
        WRITE(IUO1,70) NFIN
        DO JPLAN=1,NFIN
          IF(JPLAN.LE.NRELA) THEN
            X1=1.
            X2=0.
            PCADS=PCRELA(JPLAN)
          ELSEIF((JPLAN.GT.NRELA).AND.(JPLAN.LE.L)) THEN
            X1=0.
            X2=0.
          ELSE
            X1=0.
            X2=1.
            PCSUBS=PCREL(JPLAN-L)
          ENDIF
          REL(JPLAN)=0.
          IF(JPLAN.GT.NREL+L) GO TO 20
          IF(JPLAN.EQ.NFIN) GO TO 20
          DPLAN=VALFIN(JPLAN)-VALFIN(JPLAN+1)
          REL(JPLAN)=DPLAN*(X1*PCADS+X2*PCSUBS)/100.
  20      DREL=VALFIN(JPLAN)+REL(JPLAN)
          WRITE(IUO1,30) JPLAN,VALFIN(JPLAN),DREL
        ENDDO
C
        NBR=0
        DO JTYP=1,NAT
          NBAT=NTYP(JTYP)
          DO NUM=1,NBAT
            NBR=NBR+1
            DO JPLAN=1,NFIN
              DIF=ABS(COORD(3,NBR)-VALFIN(JPLAN))
              IF(DIF.LT.SMALL) THEN
                COORD(3,NBR)=COORD(3,NBR)+REL(JPLAN)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
        DO JPLAN=1,NFIN
          VALFIN(JPLAN)=VALFIN(JPLAN)+REL(JPLAN)
        ENDDO
C
      ELSEIF(IREL.GE.2) THEN
C
        IP=0
        LSP(2)=0
        OMEGA=OMEGA1
  97    XN1=1.
        XN2=0.
        IP=IP+1
        CALL ORDRE(NINI,VALINI,NFIN,VALFIN)
        ZP=VALFIN(IP)
        CALL RZB110(OMEGA,DY1,DY2,DZ1,DZ2)
        DZA(IP)=DZ1
        DZB(IP)=DZ2
        DYA(IP)=DY1
        DYB(IP)=DY2
        IF(ABS(OMEGA).LT.SMALL) THEN
          LSP(IP)=1
        ELSE
          LSP(IP)=2
        ENDIF
        IF(LSP(IP).EQ.1) GOTO 95
        NBR=0
C
        DO JTYP=1,NAT-NATA
          NBAT=NTYP(JTYP)
          XN1=XN1+1.-FLOAT(JTYP)
          XN2=XN2-1.+FLOAT(JTYP)
          DO JNUM=1,NBAT
            NBR=NBR+1
            ZAT=COORD(3,NBR)-ZP
            IF(ABS(ZAT).LT.SMALL) THEN
              YAT=COORD(2,NBR)
              COORD(2,NBR)=YAT-XN1*DYA(IP)-XN2*DYB(IP)
              COORD(3,NBR)=ZAT+ZP+XN1*DZA(IP)+XN2*DZB(IP)
            ENDIF
          ENDDO
        ENDDO
C
  95    OMEGA=OMEGA2
        IF((IREL.EQ.3).AND.(IP.EQ.1)) GOTO 97
        LS=0
        DO I=1,IP
          LS=LS+LSP(I)
        ENDDO
        NBZ1=NFIN+LS-IP
        DO K=1,IP
          IF(LSP(K).EQ.2) THEN
            IF((K.EQ.2).AND.(LS.EQ.3)) THEN
              KN=K-1
            ELSE
              KN=K
            ENDIF
            VALINI(NBZ1-KN+1)=VALFIN(L+K)+DZB(K)
            REL(NBZ1-KN+1)=DZB(K)
          ELSE
            VALINI(NBZ1-K+1)=VALFIN(L+K)
            REL(NBZ1-K+1)=0.
          ENDIF
        ENDDO
C
        IL=0
        IR=0
        DO J=1,NFIN
          IS=0
          IF(J.LE.NRELA) THEN
            X1=1.
            X2=0.
            X3=0.
            PCADS=PCRELA(J)
            IS=1
          ELSEIF((J.GT.NRELA).AND.(J.LE.L)) THEN
            X1=0.
            X2=0.
            X3=0.
          ELSEIF((J.GT.L).AND.(J.LE.(L+IP))) THEN
            IR=IR+1
            IF(LSP(IR).EQ.1) THEN
              IF((IR.EQ.1).AND.(LSP(2).EQ.2)) GOTO 31
              X1=0.
              X2=1.
              X3=0.
              LT=MAX0(LSP(1),LSP(2))-1
              PCSUBS=PCREL(J-L-LT)
              IL=1
              IS=1
  31          CONTINUE
            ELSE
              X1=0.
              X2=0.
              X3=1.
            ENDIF
          ELSEIF((J.GT.(L+IP)).AND.(J.LE.(L+IP+NREL))) THEN
            X1=0.
            X2=1.
            X3=0.
            LT=MAX0(LSP(1),LSP(2))+IP-1
            PCSUBS=PCREL(J-L-LT+IL+1)
            IS=1
          ELSE
            X1=0.
            X2=0.
            X3=0.
          ENDIF
          DPLAN=VALFIN(J)-VALFIN(J+1)
          REL(J)=X3*DZA(IR)+DPLAN*(X1*PCADS+X2*PCSUBS)/100.
          VALINI(J)=VALFIN(J)+REL(J)
          IF(IS.EQ.1) THEN
            NBR=0
            DO JTYP=1,NAT
              NBAT=NTYP(JTYP)
              DO NUM=1,NBAT
                NBR=NBR+1
                DIF=ABS(COORD(3,NBR)-VALFIN(J))
                IF(DIF.LT.SMALL) THEN
                  COORD(3,NBR)=VALINI(J)
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
C
        CALL ORDRE(NBZ1,VALINI,NFIN,VALFIN2)
        WRITE(IUO1,65) NFIN
        KZ(1)=0
        KZ(2)=LSP(1)
        KZ(3)=MAX0(LSP(1),LSP(2))
        DO KK=4,NFIN
          KZ(KK)=LS
        ENDDO
        DO JPLAN=1,NFIN
          IF(JPLAN.LE.L) THEN
            WRITE(IUO1,55) JPLAN,VALFIN(JPLAN),VALFIN2(JPLAN)
            VALINI(JPLAN)=VALFIN(JPLAN)
          ELSEIF((JPLAN.GT.L).AND.(JPLAN.LE.(L+LS))) THEN
            K=KZ(JPLAN-L) - INT((JPLAN-L)/2)
            IPLAN=JPLAN-K
            WRITE(IUO1,55) JPLAN,VALFIN(IPLAN),VALFIN2(JPLAN)
            VALINI(JPLAN)=VALFIN(IPLAN)
          ELSEIF(JPLAN.GT.(L+LS)) THEN
            IPLAN=JPLAN-LS+IP
            WRITE(IUO1,55) JPLAN,VALFIN(IPLAN),VALFIN2(JPLAN)
            VALINI(JPLAN)=VALFIN(IPLAN)
          ENDIF
        ENDDO
      ENDIF
C
  30  FORMAT(/,26X,'THE Z POSITION OF PLANE ',I3,' IS : ',F6.3,
     *      ' BEFORE RELAXATION AND : ',F6.3,' AFTER')
  55  FORMAT(/,26X,'THE Z  POSITION OF PLANE ',I3,' IS : ',F6.3,
     *        ' BEFORE RELAXATION AND : ',F6.3,' AFTER')
  65  FORMAT(//,44X,'THE SUMMATION IS PERFORMED OVER ',I2,' PLANES : ')
  70  FORMAT(//,44X,'THE SUMMATION IS PERFORMED OVER ',I2,' PLANES : ')
C
      RETURN
C
      END
