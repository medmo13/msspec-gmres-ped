C
C=======================================================================
C
      SUBROUTINE RZB110(OMEGA,DY1,DY2,DZ1,DZ2)
C
      A1=COS(OMEGA)
      ALPHA=SIN(OMEGA)
      BETA=A1-3.
      GAMMA=SQRT(3.)*(5./3.-A1)
      DELTA=SQRT(SQRT(3.)*(1./3.+A1)/GAMMA)
      CSA=SQRT(3.)*(-BETA-ALPHA*DELTA)/6.
      SNA=SQRT(1.-CSA*CSA)
      CSB=-SQRT(3.)*BETA/3. -CSA
      SNB=-SQRT(3.)*ALPHA/3. +SNA
      DY1=(SQRT(3.)*CSB-1.)/4.
      DY2=(1.-SQRT(3.)*CSA)/4.
      DZ1=(SQRT(3.)*SNB-SQRT(2.))/4.
      DZ2=(SQRT(3.)*SNA-SQRT(2.))/4.
C
      RETURN
C
      END
