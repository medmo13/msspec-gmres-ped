C
C=======================================================================
C
      SUBROUTINE POLLEG(NC,X,PL)
C
C  This routine computes the Legendre polynomials up to order NC-1
C
      DIMENSION PL(0:100)
C
      PL(0)=1.
      PL(1)=X
      DO 10 L=2,NC-1
        L1=L-1
        L2=L-2
        L3=2*L-1
        PL(L)=(X*FLOAT(L3)*PL(L1)-FLOAT(L1)*PL(L2))/FLOAT(L)
  10  CONTINUE
C
      RETURN
C
      END
