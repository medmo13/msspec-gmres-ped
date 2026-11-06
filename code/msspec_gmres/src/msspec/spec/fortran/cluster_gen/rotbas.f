C
C=======================================================================
C
      SUBROUTINE ROTBAS(ROT)
C
C  This routine calculates the basis vectors related to a surface
C         characterized by its Miller indices (IH,IK,II,IL)
C
      USE MILLER_MOD
      USE OUTUNITS_MOD
      USE RESEAU_MOD
      USE VECSYS_MOD , A1 => ASYS, A2 => BSYS, A3 => CSYS
C
      DIMENSION ROT(3,3),VECT(3,3),A1STAR(3),A2STAR(3),A3STAR(3),B1(3)
      DIMENSION VECT1(3),XNORM(3),CHBASE(3,3),VECT2(3,3)
C
C
C
      DATA PI /3.141593/
C
      IF((NCRIST.EQ.8).AND.(IVG0.GE.1)) GOTO 7
      XH=FLOAT(IH)
      XK=FLOAT(IK)
      XI=FLOAT(II)
      XL=FLOAT(IL)
      XI1=-XH-XK
      II1=INT(XI1)
      IF((NCRIST.EQ.7).AND.(XI.NE.XI1)) WRITE(IUO1,5) IH,IK,II1,IL
   5  FORMAT(5X,'THE SURFACE INDICES ARE NOT CORRECT,',/,5X, 'FOR THE RE
     &ST OF THE CALCULATION, THEY ARE TAKEN AS ','(',I2,1X,I2,1X,I2,1X,I
     &2,')')
      CPR=1.
      CALL PRVECT(A2,A3,B1,CPR)
      OMEGA=PRSCAL(A1,B1)/(2.*PI)
      CALL PRVECT(A2,A3,A1STAR,OMEGA)
      CALL PRVECT(A3,A1,A2STAR,OMEGA)
      CALL PRVECT(A1,A2,A3STAR,OMEGA)
      DO 10 I=1,3
        VECT1(I)=XH*A1STAR(I)+XK*A2STAR(I)+XL*A3STAR(I)
  10  CONTINUE
      DO 15 I=1,3
        ROT(I,3)=VECT1(I)/SQRT(PRSCAL(VECT1,VECT1))
  15  CONTINUE
      DO 20 I=1,3
        CHBASE(I,1)=A1(I)
        CHBASE(I,2)=A2(I)
        CHBASE(I,3)=A3(I)
        DO 25 J=1,3
          VECT(I,J)=0.
  25    CONTINUE
  20  CONTINUE
      XHKL=XH*XK*XL
      XHK=XH*XK
      XHL=XH*XL
      XKL=XK*XL
      IF(XHKL.NE.0.) THEN
        VECT(1,1)=-1./XH
        VECT(2,1)=1./XK
        VECT(1,2)=-1./XH
        VECT(3,2)=1./XL
        VECT(2,3)=-1./XK
        VECT(3,3)=1./XL
      ELSEIF(XHK.NE.0.) THEN
        VECT(1,1)=-1./XH
        VECT(2,1)=1./XK
      ELSEIF(XHL.NE.0.) THEN
        VECT(1,2)=-1./XH
        VECT(3,2)=1./XL
      ELSEIF(XKL.NE.0.) THEN
        VECT(2,3)=-1./XK
        VECT(3,3)=1./XL
      ELSEIF(XH.NE.0.) THEN
        VECT(2,2)=1./XH
      ELSEIF(XK.NE.0.) THEN
        VECT(3,3)=1./XK
      ELSEIF(XL.NE.0.) THEN
        VECT(1,1)=1./XL
      ENDIF
      CALL MULMAT(CHBASE,3,3,VECT,3,3,VECT2)
      DO 35 I=1,3
        XNORM(I)=SQRT(VECT2(1,I)**2+VECT2(2,I)**2+VECT2(3,I)**2)
  35  CONTINUE
      XMIN=AMIN1(XNORM(1),XNORM(2),XNORM(3))
      XMAX=AMAX1(XNORM(1),XNORM(2),XNORM(3))
      DO 40 I=1,3
        IF(XHKL.NE.0.) THEN
          IF(ABS(XMIN-XNORM(I)).LT.0.0001) THEN
            DO 45 J=1,3
              ROT(J,1)=VECT2(J,I)/XNORM(I)
  45        CONTINUE
          ENDIF
        ELSE
          IF(ABS(XMAX-XNORM(I)).LT.0.0001) THEN
            DO 50 J=1,3
              ROT(J,1)=VECT2(J,I)/XNORM(I)
  50        CONTINUE
          ENDIF
        ENDIF
  40  CONTINUE
      ROT(1,2)=ROT(2,3)*ROT(3,1)-ROT(3,3)*ROT(2,1)
      ROT(2,2)=ROT(3,3)*ROT(1,1)-ROT(3,1)*ROT(1,3)
      ROT(3,2)=ROT(1,3)*ROT(2,1)-ROT(2,3)*ROT(1,1)
      IF(NCRIST.EQ.7) THEN
        WRITE(IUO1,85) IH,IK,II1,IL
      ELSE
        WRITE(IUO1,80) IH,IK,IL
      ENDIF
      WRITE(IUO1,65) ROT(1,1),ROT(2,1),ROT(3,1)
      WRITE(IUO1,70) ROT(1,2),ROT(2,2),ROT(3,2)
      WRITE(IUO1,75) ROT(1,3),ROT(2,3),ROT(3,3)
      GOTO 37
   7  DO 17 I=1,3
        DO 27 J=1,3
          ROT(I,J)=0.
          IF(I.EQ.J) ROT(I,J)=1.
  27    CONTINUE
  17  CONTINUE
      IF(IVG0.EQ.1) WRITE(IUO1,48)
      IF(IVG0.EQ.2) WRITE(IUO1,47)
  47  FORMAT(//,25X,'LINEAR CHAIN STUDY ')
  48  FORMAT(//,35X,'PLANE STUDY')
  65  FORMAT(26X,'ISURF = (',F6.3,',',F6.3,',',F6.3,')')
  70  FORMAT(26X,'JSURF = (',F6.3,',',F6.3,',',F6.3,')')
  75  FORMAT(26X,'KSURF = (',F6.3,',',F6.3,',',F6.3,')')
  80  FORMAT(//,18X,'BASIS VECTORS FOR THE SURFACE (',I2,1X,I2,1X,
     *I2,') :',/)
  85  FORMAT(//,18X,'BASIS VECTORS FOR THE SURFACE (',I2,1X,I2,1X,
     *I2,1X,I2,') :',/)
C
  37  RETURN
C
      END
