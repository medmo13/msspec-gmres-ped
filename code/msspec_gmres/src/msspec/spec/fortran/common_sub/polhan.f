C
C=============================================================================
C
      SUBROUTINE POLHAN(ISPHER,NO,NC,RHO,HLM)
C
C  This routine calculates a function HLM(L,M), related to the the Hankel
C      polynomials and their derivatives with respect to z=1/ikr,
C      necessary for the Rehr-Albers expansion of the propagator.
C
      USE DIM_MOD
C
      COMPLEX HLM(0:NO_ST_M,0:NL_M-1),RHO,Z,ONEC
C
      ONEC=(1.,0.)
C
      IF(ISPHER.GE.1) THEN
        Z=(0.,-1.)/RHO
C
C  Case M = 0
C
        HLM(0,0)=ONEC
        HLM(0,1)=ONEC-Z
        DO L=2,NC
          HLM(0,L)=HLM(0,L-2)-FLOAT(L+L-1)*Z*HLM(0,L-1)
        ENDDO
C
C  Case M > 0
C
        IF(NO.GE.1) THEN
          DO M=1,NO
            HLM(M,M)=-Z*HLM(M-1,M-1)*FLOAT(M+M-1)
            HLM(M,M+1)=HLM(M,M)*FLOAT(M+M+1)*(ONEC-Z*FLOAT(M+1))
            DO L=M+2,NC
              HLM(M,L)=HLM(M,L-2)-FLOAT(L+L-1)*Z*(HLM(M,L-1)+HLM(M-1,L-1
     &))
            ENDDO
          ENDDO
        ENDIF
      ELSE
        DO M=0,NO
          DO L=M,NC
            HLM(M,L)=ONEC
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
C
      END
