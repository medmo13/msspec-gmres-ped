C
C======================================================================
C
      SUBROUTINE COEFPQ(NAT,NGR)
C
C  This subroutine computes the P(n,m) and Q(n) coefficients 
C   involved in the correlation expansion formulation
C
C  Reference  : equations (2.15) and (2.16) of 
C               H. Zhao, D. Sebilleau and Z. Wu, 
C               J. Phys.: Condens. Matter 20, 275241 (2008)
C
C                                                    H.-F. Zhao 2007
C
      USE DIM_MOD
      USE Q_ARRAY_MOD
C
      INTEGER NAT,NGR
C
      REAL CMN(NGR_M,NGR_M),P(NGR_M,NGR_M)
C
C
      IF(NGR.GT.NAT) THEN
        WRITE(6,*) 'NGR is larger than NAT, which is wrong'
        STOP
      ENDIF
C
      CALL CMNGR(NAT,NGR,CMN)
C
      DO N=1,NGR
        P(N,N)=1.
        Q(N)=P(N,N)
        DO M=N+1,NGR
          P(N,M)=0.
          DO I=N,M-1
            P(N,M)=P(N,M)-P(N,I)*CMN(I,M)
          ENDDO
          Q(N)=Q(N)+P(N,M)
C
        ENDDO
C
      ENDDO
C
      RETURN
C
      END
