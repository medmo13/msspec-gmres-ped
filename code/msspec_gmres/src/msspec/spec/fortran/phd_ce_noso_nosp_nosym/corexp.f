C
C======================================================================
C
      SUBROUTINE COREXP_SAVM(JE,IGR,NGR,NLM,ITYPE,IGS,TAU)
C
C  This subroutine call the correlation matrices calculations
C    for a given order IGR
C
C                                         H.-F. Zhao : 2007
C
      USE DIM_MOD
      USE COOR_MOD
      USE Q_ARRAY_MOD
      USE TRANS_MOD
C
      INTEGER NLM(NGR_M),ITYPE(NGR_M),IGS(NGR_M)
C
      REAL QI
C
      COMPLEX*16 TAU(LINMAX,LINFMAX,NATCLU_M)
C
C
      DO ITYP=1,N_PROT
        NBTYP=NATYP(ITYP)
        NLM(IGR)=LMAX(ITYP,JE)
        ITYPE(IGR)=ITYP
        DO NUM=1,NBTYP
          IGS(IGR)=NCORR(NUM,ITYP)
C
          IF(IGS(IGR).GT.IGS(IGR-1)) THEN
            QI=Q(IGR)
            CALL MPIS(IGR,NLM,ITYPE,IGS,JE,QI,TAU)
C
            IGR=IGR+1
            IF(IGR.LE.NGR) THEN
              CALL COREXP_SAVM1(JE,IGR,NGR,NLM,ITYPE,IGS,TAU)
            ENDIF
            IGR=IGR-1
C
          ENDIF
C
        ENDDO
      ENDDO
C
      RETURN
C
      END
