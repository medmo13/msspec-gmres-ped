C
C======================================================================
C
      SUBROUTINE COREXP_SAVM1(JE,IGR,NGR,NLM,ITYPE,IGS,TAU)
C
C  This subroutine allows a recursive use of COREXP_SAVM
C
C                                          H.-F. Zhao : 2007
C
      USE DIM_MOD
C
      INTEGER NLM(NGR_M),ITYPE(NGR_M),IGS(NGR_M)
      COMPLEX*16 TAU(LINMAX,LINFMAX,NATCLU_M)
C
      CALL COREXP_SAVM(JE,IGR,NGR,NLM,ITYPE,IGS,TAU)
C
      RETURN
C
      END
