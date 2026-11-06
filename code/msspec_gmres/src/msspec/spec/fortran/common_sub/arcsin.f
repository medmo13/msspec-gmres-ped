C
C=======================================================================
C
      SUBROUTINE ARCSIN(U,CST,RANGLE)
C
C   For a given complex number U, this subroutine calculates its phase
C   Warning : it is implicitely assumed that U = sin(theta) exp(i*phi)
C   with theta > or = to 0 which is always the case when theta is obtained
C   from the coordinates of a given vector r by the ACOS intrinsic function.
C
C   When sin(theta) = 0, then phi = 0 if cos(theta) = 1 and pi if
C   cos(theta) = -1. Cos(theta) is the variable CST.
C
      COMPLEX U,CANGLE
C
      IF(CABS(U).LT.0.0001) THEN
        IF(CST.GT.0.) THEN
          RANGLE=0.
        ELSEIF(CST.LT.0.) THEN
          RANGLE=3.141593
        ENDIF
      ELSE
        CANGLE=(0.,-1.)*CLOG(U/CABS(U))
        RANGLE=REAL(CANGLE)
      ENDIF
      RETURN
C
      END
