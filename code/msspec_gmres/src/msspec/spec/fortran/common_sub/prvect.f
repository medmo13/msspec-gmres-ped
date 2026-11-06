C
C=======================================================================
C
      SUBROUTINE PRVECT(A1,A2,A3,C)
C
C  This function computes the vector product of the two vectors A1 and A2.
C            The result is A3; C is a scaling factor
C
      DIMENSION A1(3),A2(3),A3(3)
C
      A3(1)=(A1(2)*A2(3)-A1(3)*A2(2))/C
      A3(2)=(A1(3)*A2(1)-A1(1)*A2(3))/C
      A3(3)=(A1(1)*A2(2)-A1(2)*A2(1))/C
C
      RETURN
C
      END
