C
C=======================================================================
C
      SUBROUTINE LOCATE(XX,N,X,J)
C
C
C            This subroutine is taken from the book :
C          "Numerical Recipes : The Art of Scientific
C           Computing" par W.H. PRESS, B.P. FLANNERY,
C               S.A. TEUKOLSKY et W.T. VETTERLING
C               (Cambridge University Press 1992)
C
C    It performs a search in an ordered table using a bisection method.
C    Given a monotonic array XX(1:N) and a value X, it returns J such
C                  that X is between XX(J) and XX(J+1).
C
      INTEGER J,N
      INTEGER JL,JM,JU
C
      REAL X,XX(N)
C
      JL=0
      JU=N+1
  10  IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GOTO 10
      ENDIF
      J=JL
C
      RETURN
C
      END
