C
C=======================================================================
C
      FUNCTION UJ_SQ(JTYP)
C
C  This routine evaluates the mean square displacements UJ_SQ,
C    first along une direction (x, y or z): UJ2 within the Debye model,
C              using the Debye function formulation
C
C  X1 is the Debye function phi_1
C  UJ_SQ is given in unit of the square of the lattice parameter A0
C  Temperatures are expressed in Kelvin
C
C  The coefficient COEF equals:
C
C           3 hbar^{2} N_A 10^{3} / (4 k_B)
C
C    where N_A is the Avogadro number, k_B is Boltzmann's constant
C    and 10^3 arises from the fact that the atomic mass is
C    expressed in grams
C
C  Then UJ_SQ is obtained as UJ_SQ = (2 + RSJ) UJJ for surface atoms
C                            UJ_SQ = 3 UJJ for bulk atoms
C
C
C  For empty spheres, two possibilities are provided. By construction,
C    they are very light (their mass is taken as 1/1836 of the mass
C    of a H atom) and so they will vibrate a lot (IDCM = 1). When
C    setting IDCM = 2, their mean square displacement is set to a
C    tiny value so that they hardly vibrate (frozen empty spheres)
C
C                                        Last modified : 31 Jan 2017
C
      USE DIM_MOD
C
      USE DEBWAL_MOD , T => TEMP
      USE MASSAT_MOD , XM => XMT
      USE RESEAU_MOD , N1 => NCRIST, N2 => NCENTR, N3 => IBAS, N4 => NAT
     &, A0 => A, R1 => BSURA, R2 => CSURA, UN => UNIT
      USE VIBRAT_MOD
C
      REAL MJ
C
C
C
      DATA COEF   /36.381551/                            ! 3 hbar^{2} / (4 k_B) for MJ in grams
      DATA RZ2    /1.644934/                             ! Pi^2 / 6
      DATA LITTLE /0.01/                                 ! lowest temperature for calculation of phi_1
C
      N_MAX=20
C
C  Computation of the 1D mean square displacement UJ2
C
      A=TD/T
      MJ=XM(JTYP)
      C=COEF/(MJ*TD)
C
      X1=0.
      IF(T.GT.LITTLE) THEN
        DO N=1,N_MAX
          Z=FLOAT(N)
          X1=X1+EXP(-Z*A)*(A+1./Z)/Z
        ENDDO
      ENDIF
C
      P1=1.+4.*(RZ2-X1)/(A*A)
      UJJ=C*P1/(A0*A0)
C
C  3D mean square displacement UJ_SQ
C
      IF(IDCM.EQ.1) THEN
        UJ_SQ=(3.+FLOAT(I_FREE(JTYP))*(RSJ-1.))*UJJ
      ELSEIF(IDCM.EQ.2) THEN
        UJ_SQ=1.0E-20
      ENDIF
C
      RETURN
C
      END
