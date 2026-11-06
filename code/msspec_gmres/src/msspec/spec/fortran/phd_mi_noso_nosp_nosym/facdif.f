C
C=======================================================================
C
      SUBROUTINE FACDIF(COSTH,JAT,JE,FTHETA)
C
C  This routine computes the plane wave scattering factor
C
      USE DIM_MOD
C
      USE TRANS_MOD
C
      DIMENSION PL(0:100)
C
      COMPLEX FTHETA
C
      FTHETA=(0.,0.)
      NL=LMAX(JAT,JE)+1
      CALL POLLEG(NL,COSTH,PL)
      DO 20 L=0,NL-1
        FTHETA=FTHETA+(2*L+1)*TL(L,1,JAT,JE)*PL(L)
  20  CONTINUE
      FTHETA=FTHETA/VK(JE)
C
      RETURN
C
      END
