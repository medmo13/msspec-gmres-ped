C
C=======================================================================
C
      SUBROUTINE LPM(E,XLPM,*)
C
C  This routine generates the electron mean free path
C
C       ILPM=-1: XLPM is set to 1.E+30
C       ILPM=0 : XLPM is the value given in the input data file
C       ILPM=1 : XLPM computed from Tokutaka et al, Surf. Sci. 149,349 (1985)
C       ILPM=2 : XLPM computed from the Seah and Dench expression
C
C                                          Last modified : 15 Sep 2009
C
      USE LPMOY_MOD , NZ => NZA, XMAT => XMTA, RHO => RHOTA
      USE OUTUNITS_MOD
      USE TESTS_MOD
      USE VALIN_MOD
C
      E=E+VINT
C
      IF(ILPM.EQ.-1) THEN
         XLPM=1.E+30
      ELSEIF(ILPM.EQ.0) THEN
        XLPM=XLPM0
      ELSEIF(ILPM.EQ.1) THEN
        Q=FLOAT(NZ)*RHO/XMAT
        CSTE1=ALOG(Q/4.50)/(ALOG(7.74/4.50))
        CSTE2=ALOG(Q/3.32)/(ALOG(7.74/3.32))
        CSTE3=ALOG(Q/3.32)/(ALOG(4.50/3.32))
        A1=0.7271+0.2595*ALOG(E)
        A2=-3.2563+0.9395*ALOG(E)
        A3=-2.5716+0.8226*ALOG(E)
        IF(E.GE.350.) GO TO 10
        XLN=CSTE1*(0.0107-0.0083*ALOG(E))+A1
        GO TO 20
  10    IF((NZ.GE.24).AND.(NZ.LE.74)) GO TO 30
        XLN=CSTE2*(1.6551-0.2890*ALOG(E))+A2
        GO TO 20
  30    IF(NZ.GE.42) GO TO 40
        XLN=CSTE3*(0.6847-0.1169*ALOG(E))+A2
        GO TO 20
  40    XLN=CSTE1*(0.9704-0.1721*ALOG(E))+A3
  20    XLPM=EXP(XLN)
      ELSEIF(ILPM.EQ.2) THEN
        XLPM=1430./(E**2)+0.54*SQRT(E)
      ELSE
        RETURN 1
      ENDIF
C
      E=E-VINT
      IF(IPRINT.GT.0) WRITE(IUO1,80) E,XLPM
C
  80  FORMAT(/////,2X,'=========  E = ',F7.2,' eV',5X,'MEAN',' FREE PATH
     & = ',F6.3,' ANGSTROEMS  ','=========')
C
      RETURN
C
      END
