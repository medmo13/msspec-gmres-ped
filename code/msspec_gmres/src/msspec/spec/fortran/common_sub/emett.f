C
C=======================================================================
C
      SUBROUTINE EMETT(JEM,IEMET,Z,COORD,NATYP,EMET,NM,JNEM,*)
C
C  This routine looks for the position of an absorber of type IEMET(JEM)
C    situated in the plane at Z. The result is stored in EMET(3)
C
      USE DIM_MOD
C
      DIMENSION IEMET(NEMET_M)
      DIMENSION EMET(3),DIST(NATCLU_M),COORD(3,NATCLU_M),NATYP(NATM)
C
      KEMET=0
      JNT=0
      IEM=IEMET(JEM)
      IF(IEM.GT.1) THEN
        DO JTP=1,IEM-1
          JNT=JNT+NATYP(JTP)
        ENDDO
      ENDIF
      NB=NATYP(IEM)
      XMIN=1000000.
C
      DO J=1,NB
        JN=J+JNT
        DELTAZ=ABS(COORD(3,JN)-Z)
        IF(DELTAZ.LT.0.0001) THEN
          XX=COORD(1,JN)
          XY=COORD(2,JN)
          XZ=COORD(3,JN)
          DIST(J)=SQRT(XX*XX+XY*XY+XZ*XZ)
          IF(DIST(J).LT.XMIN) THEN
            XMIN=DIST(J)
            NM=IEM
            JNEM=J
            DO I=1,3
              EMET(I)=COORD(I,JN)
            ENDDO
          ENDIF
          KEMET=KEMET+1
        ENDIF
      ENDDO
C
      IF(KEMET.EQ.0) THEN
        NM=IEM
        RETURN 1
      ENDIF
C
      RETURN
C
      END
