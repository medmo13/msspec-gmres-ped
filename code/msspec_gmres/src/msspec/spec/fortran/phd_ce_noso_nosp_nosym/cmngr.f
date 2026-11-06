C
C======================================================================
C
      SUBROUTINE CMNGR(NAT,NGR,CMN)
C
C    input   : NAT,NGR
C    output  : CMN
C
C    This subroutine calculate C(NAT-N,M-N) where,
C                    1<=M<=NGR<=NAT,1<=N<=M
C               C(NAT-N,M-N) is stored as CMN(N,M)
C
C                                              H.-F. Zhao 2007
C
      USE DIM_MOD
C
      INTEGER NAT,NGR
C
      REAL CMN(NGR_M,NGR_M)
C
      IF(NGR.GT.NAT) THEN
        WRITE(6,*) 'NGR is larger than NAT, which is wrong'
        STOP
      ENDIF
C
      DO M=1,NGR
        DO N=1,NGR
          CMN(N,M)=0.
        ENDDO
        CMN(M,M)=1.
      ENDDO
C
      DO M=1,NGR
        DO N=M-1,1,-1
          CMN(N,M)=CMN(N+1,M)*FLOAT(NAT-N)/FLOAT(M-N)
        ENDDO
      ENDDO
C
      RETURN
C
      END
