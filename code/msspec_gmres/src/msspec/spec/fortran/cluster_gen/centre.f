C
C=======================================================================
C
      SUBROUTINE CENTRE(VECBAS)
C
C  This routine modifies the Bravais lattice basis vectors according to
C                 the way the lattice is centered
C
      USE RESEAU_MOD
C
      DIMENSION VECBAS(9),V1(9)
C
C
C
      DO I=1,9
        V1(I)=VECBAS(I)
      ENDDO
      N1=NCENTR-1
      GO TO (2,3,4,5,6,7) N1
C
   2  DO I=1,3
        VECBAS(I)=-0.5*V1(I)+0.5*V1(I+3)+0.5*V1(I+6)
        VECBAS(I+3)=0.5*V1(I)-0.5*V1(I+3)+0.5*V1(I+6)
        VECBAS(I+6)=0.5*V1(I)+0.5*V1(I+3)-0.5*V1(I+6)
      ENDDO
      GO TO 8
C
   3  DO I=1,3
        VECBAS(I)=0.5*(V1(I+3)+V1(I+6))
        VECBAS(I+3)=0.5*(V1(I)+V1(I+6))
        VECBAS(I+6)=0.5*(V1(I)+V1(I+3))
      ENDDO
      GO TO 8
C
 4    DO I=1,3
        VECBAS(I)=(2./3.)*V1(I)+(1./3.)*V1(I+3)+(1./3.)*V1(I+6)
        VECBAS(I+3)=(-1./3.)*V1(I)+(1./3.)*V1(I+3)+(1./3.)*V1(I+6)
        VECBAS(I+6)=(-1./3.)*V1(I)-(2./3.)*V1(I+3)+(1./3.)*V1(I+6)
      ENDDO
      DO I=1,3
        VECBAS(3*I)=VECBAS(3*I)*SQRT(3./(1.-CSURA*CSURA))
      ENDDO
      GO TO 8
C
   5  DO I=1,3
        VECBAS(I+6)=0.5*(V1(I+3)+V1(I+6))
      ENDDO
      GO TO 8
C
   6  DO I=1,3
        VECBAS(I+6)=0.5*(V1(I)+V1(I+6))
      ENDDO
      GO TO 8
C
   7  DO I=1,3
        VECBAS(I+3)=0.5*(V1(I)+V1(I+3))
      ENDDO
C
   8  RETURN
C
      END
