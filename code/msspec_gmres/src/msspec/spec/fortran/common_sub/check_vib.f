C
C=======================================================================
C
      SUBROUTINE CHECK_VIB(NAT2)
C
C  This subroutines checks the geometrical environment of each atom
C     to identify those which can move "freely" in one direction, in
C     order to see whether the mean square displacement in this
C     direction is of bulk type or surface type
C
C  An atom is considered to move freely in one direction if no other
C     atom is present in the tetragonal cell of height ALENGTH * A
C     and base edge 2 * A, whose base is centered on the atom considered
C
C  Only prototypical atoms are considered as all equivalent atoms are
C     in the same geometrical environment
C
C  Surface-like atoms are then identified as having I_FREE = 1
C
C                                     Last modified : 24 Apr 2013
C
      USE DIM_MOD
C
      USE COOR_MOD , COORD => SYM_AT
      USE OUTUNITS_MOD
      USE VIBRAT_MOD
C
      INTEGER NSUR(NATP_M)
C
      DATA SMALL /0.0001/
C
      ALENGTH=4.
C
C.................... Checking the z direction ....................
C
      WRITE(IUO1,11)
      N_SUR=0
C
C  Loop on the prototypical atoms
C
      DO JTYP=1,N_PROT
C
        I_FREE(JTYP)=0
        JAT0=NCORR(1,JTYP)
        XA=COORD(1,JAT0)
        YA=COORD(2,JAT0)
        ZA=COORD(3,JAT0)
C
C  Loop on the surrounding atoms
C
        I_ACC=0
C
        DO JAT=1,NAT2
C
          IF(JAT.EQ.JAT0) GOTO 10
C
          X=COORD(1,JAT)
          Y=COORD(2,JAT)
          Z=COORD(3,JAT)
C
C  Considering only atoms with Z > ZA
C
          IF(Z.LT.(ZA+SMALL)) GOTO 10
C
C  Lateral and vertical distances between the two atoms
C
          D_LAT=(X-XA)*(X-XA)+(Y-YA)*(Y-YA)
          D_VER=(Z-ZA)*(Z-ZA)
C
          IF(D_VER.LT.(ALENGTH+SMALL)) THEN
            IF(D_LAT.LT.(1.+SMALL)) THEN
              I_ACC=I_ACC+1
            ENDIF
          ENDIF
C
          IF(I_ACC.GE.1) GOTO 10
C
  10      CONTINUE
C
        ENDDO
C
        IF(I_ACC.EQ.0) THEN
          I_FREE(JTYP)=1
          N_SUR=N_SUR+1
          NSUR(N_SUR)=JTYP
         ENDIF
C
      ENDDO
C
      WRITE(IUO1,12) (NSUR(J),J=1,N_SUR)
C
  11  FORMAT(//,18X,'SURFACE-LIKE ATOMS FOR MSD CALCULATIONS: ',/)
  12  FORMAT(20X,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5)
C
      RETURN
C
      END
