C
C=======================================================================
C
      SUBROUTINE REFRAC(VINT,EKIN,RTHETA,RTHINT)
C
C  This routine calculates the refraction of a plane wave beam induced
C     by the surface potential barrier VINT. EKIN is the kinetic energy
C     outside the crystal.
C
C                                          Last modified :  3 Dec 2008
C
      DATA PIS180,SMALL /0.017453,0.001/
C
      IF(VINT.LT.0.) VINT=ABS(VINT)
      IF(ABS(VINT).LT.SMALL) THEN
        RTHINT=RTHETA
      ELSE
        U=VINT/(EKIN+VINT)
        DTHETA=RTHETA/PIS180
        REFRA=SIN(RTHETA)*SIN(RTHETA)*(1.-U)
        RTHINT=ASIN(SQRT(REFRA))
        IF(DTHETA.LT.0.) THEN
          RTHINT=-RTHINT
        ENDIF
      ENDIF
C
      RETURN
C
      END
