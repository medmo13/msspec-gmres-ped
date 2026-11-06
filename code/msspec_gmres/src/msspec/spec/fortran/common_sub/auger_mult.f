C
C=======================================================================
C
      SUBROUTINE AUGER_MULT
C
C  This subroutine computes all the possible multiplets that are
C     contained in a given Auger transition line. It assumes that
C     the atom has closed shells only.
C
C                                        Last modified :  9 March 2006
C
      USE INIT_A_MOD , LI => LI_C, L2 => LI_I, L1 => LI_A
      USE OUTUNITS_MOD
C
      CHARACTER*1 SC(0:1),LC(0:6),JC(0:7)
      CHARACTER*3 MULTIPLET(112)
C
      DATA SC /'1','3'/
      DATA LC /'S','P','D','F','G','H','I'/
      DATA JC /'0','1','2','3','4','5','6','7'/
C
      WRITE(IUO1,10)
      N_MULT=0
      DO NS=0,1
        DO L=ABS(L1-L2),L1+L2
          DO J=ABS(L-NS),L+NS
            N_MULT=N_MULT+1
            MULTIPLET(N_MULT)=SC(NS)//LC(L)//JC(J)
            WRITE(IUO1,20) MULTIPLET(N_MULT)
          ENDDO
        ENDDO
      ENDDO
C
  10  FORMAT(///,26X,'THE POSSIBLE MULTIPLETS ARE :',/,'  ')
  20  FORMAT(58X,A3)
C
      RETURN
C
      END
