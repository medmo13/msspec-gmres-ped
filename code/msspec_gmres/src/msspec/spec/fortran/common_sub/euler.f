C
C=======================================================================
C
      SUBROUTINE EULER(RTHETA1,RPHI1,RTHETA2,RPHI2,RALPHA,RBETA,RGAMMA,I
     &ROT)
C
C  This routine calculates the Euler angles RALPHA,RBETA,RGAMMA corresponding
C       to the rotation r1(RTHETA1,RPHI1) ----> r2(RTHETA2,RPHI2)
C
C       IROT=1 : r ---> z represented by (0,RTHETA,PI-RPHI)
C       IROT=0 : r ---> z represented by (0,-RTHETA,-RPHI)
C
C
      COMPLEX U1,U2
C
      DATA PI /3.141593/
C
      IF(IROT.EQ.1) THEN
        EPS=1
      ELSE
        EPS=-1
      ENDIF
      DPHI=RPHI2-RPHI1
      A1=SIN(RTHETA1)*COS(RTHETA2)
      A2=COS(RTHETA1)*SIN(RTHETA2)
      A3=COS(RTHETA1)*COS(RTHETA2)
      A4=SIN(RTHETA1)*SIN(RTHETA2)
      U1=A1-A2*COS(DPHI)-(0.,1.)*SIN(RTHETA2)*SIN(DPHI)
      U2=A1*COS(DPHI)-A2+(0.,1.)*SIN(RTHETA1)*SIN(DPHI)
      U3=A3+A4*COS(DPHI)
      IF(U3.GT.1.) U3=1.
      IF(U3.LT.-1.) U3=-1.
      RBETA=ACOS(U3)
      IF(ABS(SIN(RBETA)).GT.0.0001) THEN
        U1=EPS*U1/SIN(RBETA)
        U2=EPS*U2/SIN(RBETA)
        CALL ARCSIN(U1,U3,RALPHA)
        CALL ARCSIN(U2,U3,RGAMMA)
      ELSE
        RALPHA=0.
        IF(ABS(U3-1.0).LT.0.0001) THEN
          RGAMMA=0.
        ELSE
          RGAMMA=PI
        ENDIF
      ENDIF
C
      RETURN
C
      END
