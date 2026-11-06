C
C=======================================================================
C
      SUBROUTINE CONV_SERIES(ACC,N_TABLE,JORD,R_SN,*)
C
C  This subroutine checks the convergence of the power method or of
C  the Rayleigh quotient method at level JORD without convergence
C  acceleration. For convergence, there must be N_TABLE equal values
C             (ACC being the accepted error) of SN
C
C  Author : D. Sebilleau
C
C                                 Last modified : 27 Feb 2013
C
C      INCLUDE 'spec.inc'
C
      USE DIM_MOD
C
      REAL ACC,COMP,R_SN(N_ORD_M)
C
      NCONV=0
      DO J=JORD,JORD-N_TABLE,-1
        COMP=ABS(R_SN(J)-R_SN(J-1))
        IF(COMP.LE.ACC) THEN
          NCONV=NCONV+1
        ENDIF
      ENDDO
C
      IF(NCONV.EQ.N_TABLE-1) RETURN 1
C
      RETURN
C
      END
