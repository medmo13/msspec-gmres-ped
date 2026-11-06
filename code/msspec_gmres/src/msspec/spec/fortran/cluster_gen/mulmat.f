C
C=======================================================================
C
      SUBROUTINE MULMAT(A1,IL1,IC1,A2,IL2,IC2,A3)
C
C  This routine performs the matrix multiplication of A1(IL1,IC1) by
C            A2(IL2,IC2) with the result stored in A3(IL1,IC2)
C
      USE OUTUNITS_MOD
      DIMENSION A1(IL1,IC1),A2(IL2,IC2),A3(IL1,IC2)
C
C
      IF(IC1.NE.IL2) THEN
        WRITE(IUO1,10)
      ELSE
        DO I=1,IL1
          DO J=1,IC2
            A3(I,J)=0.
            DO K=1,IC1
              A3(I,J)=A3(I,J)+A1(I,K)*A2(K,J)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
  10  FORMAT(5X,'THESE MATRICES CANNOT BE MULTIPLIED')
C
      RETURN
C
      END
