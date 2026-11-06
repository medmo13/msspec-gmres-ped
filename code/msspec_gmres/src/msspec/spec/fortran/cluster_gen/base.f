C
C=======================================================================
C
      SUBROUTINE BASE
C
C  This routine generates the lattice basis vectors for a given Bravais
C           lattice NCRIST centered according to NCENTR
C
      USE DIM_MOD
      USE BASES_MOD
      USE CRANGL_MOD
      USE OUTUNITS_MOD
      USE RESEAU_MOD
      USE VECSYS_MOD
C
      CHARACTER*15  BRAV(8),CENT(7)
      CHARACTER*31 RESEAU
C
C
      DIMENSION CUB(9),MNC(9),TCN(9),TRG(9),HEX(9)
C
C
C
      DATA CUB /1.,0.,0.,   0.,1.,0.,   0.,0.,1./
      DATA MNC /1.,0.,1.,   0.,1.,0.,   0.,0.,1./
      DATA TCN /1.,0.,1.,   1.,1.,1.,   0.,0.,1./
      DATA TRG /0.,1.,1.,   -0.866025,-0.5,1.,   0.866025,-0.5,1./
      DATA HEX /1.,0.,0.,   -0.5,0.866025,0.,  0.,0.,1./
      DATA PIS180 /0.017453/
      DATA BRAV /'        CUBIQUE','     TETRAGONAL',' ORTHORHOMBIQUE','
     &   MONOCLINIQUE','    TRICLINIQUE','       TRIGONAL','      HEXAGO
     &NAL','        EXTERNE'/
      DATA CENT /' ','CENTRE',' FACES CENTREES','(RHOMBOEDRIQUE)',' FACE
     & A CENTREE',' FACE B CENTREE',' FACE C CENTREE'/
C
      ALPHAR=ALPHAD*PIS180
      BETAR=BETAD*PIS180
      GAMMAR=GAMMAD*PIS180
      NAT3=NAT*3
      GO TO (1,1,1,2,3,4,5,6) NCRIST
C
   1  DO I=1,9
        VECBAS(I)=CUB(I)
      ENDDO
      IF(NCRIST.NE.1) THEN
        VECBAS(9)=CSURA
        IF(NCRIST.EQ.3) THEN
          VECBAS(5)=BSURA
        ENDIF
      ENDIF
      GO TO 6
C
   2  DO I=1,9
        VECBAS(I)=MNC(I)
      ENDDO
      VECBAS(1)=SIN(BETAR)
      VECBAS(3)=COS(BETAR)
      VECBAS(5)=BSURA
      VECBAS(9)=CSURA
      GO TO 6
C
   3  DO I=1,9
        VECBAS(I)=TCN(I)
      ENDDO
      VECBAS(1)=SIN(BETAR)
      VECBAS(3)=COS(BETAR)
      A2Y=(COS(GAMMAR)-COS(ALPHAR)*COS(BETAR))/SIN(BETAR)
      VECBAS(4)=BSURA*A2Y
      VECBAS(5)=BSURA*SQRT(SIN(ALPHAR)*SIN(ALPHAR)-A2Y*A2Y)
      VECBAS(6)=BSURA*COS(ALPHAR)
      VECBAS(9)=CSURA
      GO TO 6
C
   4  IF(((NCENTR.EQ.4).AND.(CSURA.NE.1.)).OR.(NCENTR.EQ.1)) GO TO 5
      ETA=-2.*SIN(ALPHAR/2.)/SQRT(3.)
      DZETA=SQRT(1.-ETA*ETA)
      DO I=1,3
        J=I+2*(I-1)
        J1=J+1
        J2=J+2
        VECBAS(J)=TRG(J)*ETA
        VECBAS(J1)=TRG(J1)*ETA
        VECBAS(J2)=TRG(J2)*DZETA
      ENDDO
      GO TO 6
C
   5  DO I=1,9
        VECBAS(I)=HEX(I)
      ENDDO
      VECBAS(9)=CSURA
C
   6  DO I=1,3
        ASYS(I)=VECBAS(I)
        BSYS(I)=VECBAS(I+3)
        CSYS(I)=VECBAS(I+6)
      ENDDO
      DCA=ABS(CSURA-1.)
      IF((NCRIST.EQ.6).AND.(DCA.LT.0.0001)) GO TO 8
      IF(NCRIST.EQ.8) GO TO 8
      IF(NCENTR.GT.1) THEN
        CALL CENTRE(VECBAS)
        IF(NCENTR.EQ.4) THEN
          DO I=1,9
            VECBAS(I)=VECBAS(I)*SQRT((1.-CSURA*CSURA)*3.)
          ENDDO
          DO I=1,3
            ASYS(I)=VECBAS(I)
            BSYS(I)=VECBAS(I+3)
            CSYS(I)=VECBAS(I+6)
          ENDDO
        ENDIF
      ENDIF
C
   8  RESEAU=BRAV(NCRIST)//' '//CENT(NCENTR)
      WRITE(IUO1,80) RESEAU,NAT
      WRITE(IUO1,81) (VECBAS(I),I=1,9)
      WRITE(IUO1,82)
      WRITE(IUO1,83) (ATBAS(I),I=1,NAT3)
C
  80  FORMAT(////,10X,'RESEAU CRISTALLIN DE TYPE : ',A29,/,16X,
     *       'CONTENANT',I3,' ATOMES DANS LA MAILLE ELEMENTAIRE',//)
  81  FORMAT(28X,'VECTEURS GENERATEURS :',//,26X,'A1 = (',F6.3,',',
     *F6.3,',',F6.3,')',/,26X,'A2 = (',F6.3,',',F6.3,',',F6.3,')',/,
     *26X,'A3 = (',F6.3,',',F6.3,',',F6.3,')')
  82  FORMAT(/,21X,'POSITIONS DES ATOMES DANS LA MAILLE :',/)
  83  FORMAT(29X,'(',F6.3,',',F6.3,',',F6.3,')')
C
      RETURN
C
      END
