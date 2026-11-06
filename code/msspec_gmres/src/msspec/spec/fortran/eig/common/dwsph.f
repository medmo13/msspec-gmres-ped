C
C=======================================================================
C
      SUBROUTINE DWSPH(JTYP,JE,X,TLT,ISPEED)
C
C  This routine recomputes the T-matrix elements taking into account the
C                mean square displacements.
C
C  When the argument X is tiny, no vibrations are taken into account
C
C                                      Last modified : 25 Apr 2013
C
C      INCLUDE 'spec.inc'
      USE DIM_MOD
      USE TRANS_MOD
C
      DIMENSION GNT(0:N_GAUNT)
C
      COMPLEX TLT(0:NT_M,4,NATM,NE_M),SL1,ZEROC
C
      COMPLEX*16 FFL(0:2*NL_M)
C
      DATA PI4,EPS /12.566371,1.0E-10/
C
      ZEROC=(0.,0.)
C
      IF(X.GT.EPS) THEN
C
C  Standard case: vibrations
C
        IF(ISPEED.LT.0) THEN
          NSUM_LB=ABS(ISPEED)
        ENDIF
C
        COEF=PI4*EXP(-X)
        NL2=2*LMAX(JTYP,JE)+2
        IBESP=5
        MG1=0
        MG2=0
C
        CALL BESPHE(NL2,IBESP,X,FFL)
C
        DO L=0,LMAX(JTYP,JE)
          XL=FLOAT(L+L+1)
          SL1=ZEROC
C
          DO L1=0,LMAX(JTYP,JE)
            XL1=FLOAT(L1+L1+1)
            CALL GAUNT(L,MG1,L1,MG2,GNT)
            L2MIN=ABS(L1-L)
            IF(ISPEED.GE.0) THEN
              L2MAX=L1+L
            ELSEIF(ISPEED.LT.0) THEN
              L2MAX=L2MIN+2*(NSUM_LB-1)
            ENDIF
            SL2=0.
C
            DO L2=L2MIN,L2MAX,2
              XL2=FLOAT(L2+L2+1)
              C=SQRT(XL1*XL2/(PI4*XL))
              SL2=SL2+C*GNT(L2)*REAL(DREAL(FFL(L2)))
            ENDDO
C
            SL1=SL1+SL2*TL(L1,1,JTYP,JE)
          ENDDO
C
          TLT(L,1,JTYP,JE)=COEF*SL1
C
        ENDDO
C
      ELSE
C
C  Argument X tiny: no vibrations
C
        DO L=0,LMAX(JTYP,JE)
C
          TLT(L,1,JTYP,JE)=TL(L,1,JTYP,JE)
C
        ENDDO
C
      ENDIF
C
      RETURN
C
      END
