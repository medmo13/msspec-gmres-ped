C
C=======================================================================
C
      SUBROUTINE CHBASE(NAT,ATBAS)
C
      USE VECSYS_MOD
C
      DIMENSION ATBAS(3*NAT),BASVEC(3,3),BAS1(1,3),BAS2(1,3)
C
      DO J=1,3
        BASVEC(1,J)=ASYS(J)
        BASVEC(2,J)=BSYS(J)
        BASVEC(3,J)=CSYS(J)
      ENDDO
C
      DO JAT=1,NAT
        DO J=1,3
          K=J+3*(JAT-1)
          BAS1(1,J)=ATBAS(K)
        ENDDO
        CALL MULMAT(BAS1,1,3,BASVEC,3,3,BAS2)
        DO J=1,3
          K=J+3*(JAT-1)
          ATBAS(K)=BAS2(1,J)
        ENDDO
      ENDDO
C
      RETURN
C
      END
