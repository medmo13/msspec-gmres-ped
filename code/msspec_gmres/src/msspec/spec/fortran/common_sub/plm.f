C
C=======================================================================
C
      SUBROUTINE PLM(X,PLMM,NC)
C
C  This routine computes the Legendre functions. It is a modified version
C       of that written by W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY and
C       W.T. VETTERLING in "Numerical Recipes : The Art of Scientific
C           Computing" (Cambridge University Press 1992).
C
      DIMENSION PLMM(0:100,0:100)
C
      PLMM(0,0)=1.
      PLMM(1,0)=X
      DO L=2,NC
        PLMM(L,0)=(X*(L+L-1)*PLMM(L-1,0)-(L-1)*PLMM(L-2,0))/L
      ENDDO
C
      DO M=1,NC
        PMM=1.
        FACT=1.
        SOMX2=SQRT(1.-X*X)
        FACT=1.
        DO I=1,M
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+2.
        ENDDO
        PMMP1=X*FACT*PMM
        PLMM(M,M)=PMM
        PLMM(M+1,M)=PMMP1
        IF(M.LT.NC-1) THEN
          DO L=M+2,NC
            PLL=(X*(L+L-1)*PMMP1-(L+M-1)*PMM)/(L-M)
            PMM=PMMP1
            PMMP1=PLL
            PLMM(L,M)=PLL
          ENDDO
        ENDIF
      ENDDO
C
      RETURN
C
      END
