C
C=======================================================================
C
      SUBROUTINE INTERP_POINTS(I_XN,N_APP,ALPHA,BETA,SN,XN)
C
C  This subroutine computes the interpolation
C     points used in some algorithms
C
C
C  Author : D. Sebilleau
C
C                                 Last modified : 27 Feb 2013
C
C      INCLUDE 'spec.inc'
C
      USE DIM_MOD
C
      COMPLEX*16 SN(-1:N_ORD_M,-1:N_ORD_M),XN(-1:N_ORD_M)
      COMPLEX*16 ALPHA,BETA
C
      DO N=0,N_APP
C
        IF(I_XN.EQ.1) THEN
          XN(N)=(SN(N+1,0)-SN(N,0))
        ELSEIF(I_XN.EQ.2) THEN
          XN(N)=BETA+DFLOAT(N)
        ELSEIF(I_XN.EQ.3) THEN
          XN(N)=(BETA+DFLOAT(N))*(SN(N+1,0)-SN(N,0))
        ELSEIF(I_XN.EQ.4) THEN
          XN(N)=(BETA+DFLOAT(N))**ALPHA
        ELSEIF(I_XN.EQ.5) THEN
          XN(N)=((2.D0**N)/BETA)**ALPHA
        ENDIF
C
      ENDDO
C
      RETURN
C
      END
