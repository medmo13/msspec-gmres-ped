C
C=======================================================================
C
      SUBROUTINE PLOTFD(A,LMX,ITL,NL,NAT,NE)
C
C  This routine prepares the output for a plot of the scattering factor
C
C      INCLUDE 'spec.inc'
C
      USE DIM_MOD
C
      USE APPROX_MOD
      USE FDIF_MOD
      USE INIT_L_MOD, L => LI, I2 => INITL, I3 => NNL, I4 => LF1, I5 => 
     &LF2, I10 => ISTEP_LF
      USE INIT_J_MOD
      USE OUTFILES_MOD
      USE OUTUNITS_MOD
      USE PARCAL_MOD, N3 => NPHI, N4 => NE, N5 => NTHETA, N6 => NEPS
      USE TYPCAL_MOD, I7 => IFTHET, I8 => IMOD, I9 => IPOL, I12 => I_CP,
     & I13 => I_EXT, I14 => I_TEST
      USE VALIN_MOD, U1 => THLUM, U2 => PHILUM, U3 => ELUM, N7 => NONVOL
      USE VALFIN_MOD
C
      DIMENSION LMX(NATM,NE_M)
C
      COMPLEX FSPH,VKE
C
      DATA PI,CONV/3.141593,0.512314/
C
      OPEN(UNIT=IUO3, FILE=OUTFILE3, STATUS='UNKNOWN')
      IF(ISPHER.EQ.0) THEN
         L=0
         LMAX=0
      ELSE
        LMAX=L
      ENDIF
      PHITOT=360.
      THTOT=360.*ITHETA*(1-IPHI)+180.*ITHETA*IPHI
      NPHI=(NFTHET+1)*IPHI+(1-IPHI)
      NTHT=(NFTHET+1)*ITHETA*(1-IPHI)+(NFTHET/2+1)*ITHETA*IPHI+
     *     (1-ITHETA)
      NE=NFTHET*IE + (1-IE)
      WRITE(IUO3,1) ISPHER,NL,NAT,L,NTHT,NPHI,NE,E0,EFIN
      DO 10 JT=1,NTHT
        DTHETA=THETA1+FLOAT(JT-1)*THTOT/FLOAT(MAX0(NTHT-1,1))
        RTHETA=DTHETA*PI/180.
        TEST=SIN(RTHETA)
        IF(TEST.GE.0.) THEN
          POZ=PI
          EPS=1.
        ELSE
          POZ=0.
          EPS=-1.
        ENDIF
        BETA=RTHETA*EPS
        IF(ABS(TEST).LT.0.0001) THEN
          NPHIM=1
        ELSE
          NPHIM=NPHI
        ENDIF
        DO 20 JP=1,NPHIM
          DPHI=PHI1+FLOAT(JP-1)*PHITOT/FLOAT(MAX0(NPHI-1,1))
          RPHI=DPHI*PI/180.
          GAMMA=POZ-RPHI
          DO 30 JE=1,NE
            IF(NE.EQ.1) THEN
              ECIN=E0
            ELSE
              ECIN=E0+FLOAT(JE-1)*(EFIN-E0)/FLOAT(NE-1)
            ENDIF
            IF(ITL.EQ.0) VKE=SQRT(ECIN-ABS(VINT))*CONV*A*(1.,0.)
            DO 40 JAT=1,NAT
              IF(L.GT.LMX(JAT,JE)) GOTO 90
              DO 50 M=-LMAX,LMAX
                CALL FACDIF1(VKE,R1,R2,THETA0,PHI0,BETA,
     &          GAMMA,L,M,FSPH,JAT,JE,*60)
                GOTO 70
  60            WRITE(IUO1,80)
                STOP
  70            REFTH=REAL(FSPH)
                XIMFTH=AIMAG(FSPH)
                WRITE(IUO3,5) JE,JAT,L,M,REFTH,XIMFTH,DTHETA,
     &          DPHI,ECIN
  50          CONTINUE
              GOTO 40
  90          WRITE(IUO1,100) JAT
              STOP
  40        CONTINUE
  30      CONTINUE
  20    CONTINUE
  10  CONTINUE
      CLOSE(IUO3)
   1  FORMAT(5X,I1,2X,I2,2X,I4,2X,I2,2X,I3,2X,I3,2X,I3,2X,F8.2,2X,F8.2)
   5  FORMAT(1X,I3,1X,I4,1X,I2,1X,I3,1X,F6.3,1X,F6.3,1X,F6.2,1X,F6.2,
     &1X,F8.2)
  80  FORMAT(15X,'<<<<<  WRONG VALUE OF THETA0 : THE DENOMINATOR ','IS 
     &ZERO  >>>>>')
 100  FORMAT(15X,'<<<<<  THE VALUE OF L EST IS TOO LARGE FOR ATOM',' : 
     &',I2,'  >>>>>')
C
      RETURN
C
      END
