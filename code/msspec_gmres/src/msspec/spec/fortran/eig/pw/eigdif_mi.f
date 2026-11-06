C
C=======================================================================
C
      SUBROUTINE EIGDIF_MI
C
C   This subroutine computes the eigenvalues of the
C     multiple scattering matrix to obtain its spectral radius
C
C   The eigenvalue calculation is performed using the LAPACK
C     eigenvalue routines for a general complex matrix (I_PWM = 0)
C     or using the power method  (I_PWM > 0)
C
C                                      Last modified : 26 Apr 2013
C
C      INCLUDE 'spec.inc'
      USE DIM_MOD
      USE CONVTYP_MOD
      USE COOR_MOD, NTCLU => NATCLU, NTP => NATYP
      USE DEBWAL_MOD
      USE EIGEN_MOD, NE => NE_EIG, E0 => E0_EIG, EFIN => EFIN_EIG
      USE OUTFILES_MOD
      USE OUTUNITS_MOD
      USE RESEAU_MOD
      USE TESTS_MOD
      USE TRANS_MOD
      USE VALIN_MOD, E1 => E0
C
      COMPLEX IC,ONEC
      COMPLEX TLT(0:NT_M,4,NATM,NE_M)
C
      DATA CONV /0.512314/
C
      IC=(0.,1.)
      ONEC=(1.,0.)
C
      OPEN(UNIT=IUO2, FILE=OUTFILE2, STATUS='UNKNOWN')
C
C  Check file to monitor the convergence of the calculation
C
      OPEN(UNIT=36, FILE='div/radius_conv.lis', STATUS='unknown')
C
C   Loop over the energies
C
      DO JE=1,NE
        IF(NE.GT.1) THEN
          ECIN=E0+FLOAT(JE-1)*(EFIN-E0)/FLOAT(NE-1)
        ELSEIF(NE.EQ.1) THEN
          ECIN=E0
        ENDIF
        CALL LPM(ECIN,XLPM,*6)
        XLPM1=XLPM/A
        IF(IPRINT.GT.0) WRITE(IUO1,56) A,XLPM1
        IF(ITL.EQ.0) THEN
          VK(JE)=SQRT(ECIN+VINT)*CONV*A*ONEC
          VK2(JE)=CABS(VK(JE)*VK(JE))
        ENDIF
        GAMMA=1./(2.*XLPM1)
        IF(IPOTC.EQ.0) THEN
          VK(JE)=VK(JE)+IC*GAMMA
        ENDIF
        IF(I_MFP.EQ.0) THEN
          VK(JE)=CMPLX(REAL(VK(JE)))
          VK2(JE)=CABS(VK(JE)*VK(JE))
        ENDIF
        IF(I_VIB.EQ.1) THEN
          IF(IDCM.GE.1) WRITE(IUO1,22)
          DO JAT=1,N_PROT
            IF(IDCM.EQ.0) THEN
              XK2UJ2=VK2(JE)*UJ2(JAT)
            ELSE
              XK2UJ2=VK2(JE)*UJ_SQ(JAT)
              WRITE(IUO1,23) JAT,UJ_SQ(JAT)*A*A
            ENDIF
            CALL DWSPH(JAT,JE,XK2UJ2,TLT,I_VIB)
            DO LAT=0,LMAX(JAT,JE)
              TL(LAT,1,JAT,JE)=TLT(LAT,1,JAT,JE)
            ENDDO
          ENDDO
        ENDIF
C
        IF(JE.EQ.1) WRITE(36,57)
        WRITE(36,59)
        WRITE(36,58) ECIN
        WRITE(36,59)
C
C     Eigenvalue calculation by power method
C
        CALL SPEC_RAD_POWER(JE,ECIN)
C
        WRITE(36,59)
        WRITE(36,57)
C
C     End of the loop on the energy
C
      ENDDO
      CLOSE(36)
      GOTO 7
C
   6  WRITE(IUO1,55)
C
  22  FORMAT(16X,'INTERNAL CALCULATION OF MEAN SQUARE DISPLACEMENTS',/,
     &25X,' BY DEBYE UNCORRELATED MODEL:',/)
  23  FORMAT(21X,'ATOM TYPE ',I5,'  MSD = ',F8.6,' ANG**2')
  55  FORMAT(///,12X,' <<<<<<<<<<  THIS VALUE OF ILPM IS NOT',
     &'AVAILABLE  >>>>>>>>>>')
  56  FORMAT(4X,'LATTICE PARAMETER A = ',F6.3,' ANGSTROEMS',4X,
     *'MEAN FREE PATH = ',F6.3,' * A',//)
  57  FORMAT('+++++++++++++++++++++++++++++++++++++++++++++++++++++++',
     &'++++')
  58  FORMAT(18X,'E_KIN = ',F8.2,1X,'eV')
  59  FORMAT(20X)
C
   7  RETURN
C
      END
