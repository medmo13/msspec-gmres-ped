C
C=======================================================================
C
      SUBROUTINE CHNOT(NVEC,VEC1,VEC2)
C
C  This routine linearizes the storage of a two index array
C
      DIMENSION VEC1(3*NVEC),VEC2(3,NVEC)
C
      DO J=1,NVEC
        DO I=1,3
          VEC2(I,J)=VEC1(I+3*(J-1))
        ENDDO
      ENDDO
C
      RETURN
C
      END
