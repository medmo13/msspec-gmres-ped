C
C=======================================================================
C
      SUBROUTINE INVMAT(B,BINV)
C
      USE OUTUNITS_MOD
      DIMENSION B(3,3),BINV(3,3)
C
C
      A1=B(1,1)*B(2,2)*B(3,3)
      A2=B(2,1)*B(3,2)*B(1,3)
      A3=B(3,1)*B(1,2)*B(2,3)
      A4=B(1,1)*B(3,2)*B(2,3)
      A5=B(2,1)*B(1,2)*B(3,3)
      A6=B(3,1)*B(2,2)*B(1,3)
      DET=A1+A2+A3-A4-A5-A6
C
      IF(ABS(DET).LT.0.0001) GO TO 10
C
      DO I=1,3
        DO J=1,3
          DO K=1,3
            L=(I-J)*(I-K)*(J-K)
            IF(L.NE.0) THEN
              XNUM1=B(J,J)*B(K,K)-B(J,K)*B(K,J)
              XNUM2=B(I,K)*B(K,J)-B(K,K)*B(I,J)
              BINV(I,I)=XNUM1/DET
              BINV(I,J)=XNUM2/DET
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      GO TO 50
C
  10  WRITE(IUO1,60)
C
  60  FORMAT(5X,'NON INVERTIBLE MATRIX')
C
  50  CONTINUE
C
      RETURN
C
      END
