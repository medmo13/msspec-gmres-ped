C
C=======================================================================
C
      SUBROUTINE FINDPATHS3(ND,ITYP,IATL,I_CP,R,XR,YR,ZR,RHOMI,THMI,
     &                      PHIMI,ZSURF,JPOS,PW,JE,FREF,DIJ,TAU)
C
C  This routine generates all the paths and filters them according to the
C      criteria given in the input data file (IFSPH,IFWD,IPW,ILENGTH).
C      It corresponds to the spin-independent case from a non spin-orbit
C                    resolved initial core state LI
C
C                                          Last modified : 16 May 2007
C
      USE DIM_MOD
C
      USE APPROX_MOD , ILE => ILENGTH, RLE => RLENGTH
      USE COOR_MOD
      USE DEBWAL_MOD
      USE INIT_L_MOD
      USE PATH_MOD
      USE ROT_MOD
      USE TESTPA_MOD
      USE TESTPB_MOD
      USE TRANS_MOD
      USE TLDW_MOD
      USE VARIA_MOD
C
      DIMENSION XR(NDIF_M),YR(NDIF_M),ZR(NDIF_M)
      DIMENSION JPOS(NDIF_M,3),R(NDIF_M)
C
C
C
      COMPLEX PW1,PWI,FTHETA,RHOMI,RHOIJ,RHOJK
      COMPLEX IC,COMPL1,PW(0:NDIF_M)
      COMPLEX TAU(LINMAX,LINFMAX,NATCLU_M)
      COMPLEX YLM1(0:NL_M,-NL_M:NL_M)
      COMPLEX YLM2(0:NL_M,-NL_M:NL_M),CTL,CTL2
C
      DATA XCOMP,PI4,SMALL /1.E-10,12.566371,0.0001/
C
      IC=(0.,1.)
      IEULER=1
C
      IF(IFWD.EQ.1) COSFWDI=COS(RTHFWD(ITYP))
      IF(IBWD(ITYP).EQ.1) COSBWDI=COS(RTHBWD(ITYP))
C
C    I_CP = 0 : all open paths generated
C    I_CP = 1 : only closed paths generated
C
      IF((I_CP.EQ.0).OR.(ND.NE.NDIF-1)) THEN
        N_TYP=N_PROT
      ELSE
        N_TYP=1
      ENDIF
C
      DO JTYP=1,N_TYP
       IF(IFWD.EQ.1) COSFWDJ=COS(RTHFWD(JTYP))
       IF(IBWD(JTYP).EQ.1) COSBWDJ=COS(RTHBWD(JTYP))
       ND=ND+1
C
C     I_ABS = 0 : the atom before the scatterer is not the absorber
C     I_ABS = 1 : the atom before the scatterer is the absorber
C     I_ABS = 2 : the atom after the scatterer is the absorber (XAS only)
C
       IF(ND.EQ.1) THEN
         I_ABS=1
       ELSE
         I_ABS=0
       ENDIF
C
       IF((I_CP.EQ.0).OR.(ND.NE.NDIF)) THEN
         NBTYPJ=NATYP(JTYP)
       ELSE
         NBTYPJ=1
       ENDIF
C
       DO JNUM=1,NBTYPJ
          JATL=NCORR(JNUM,JTYP)
          IF(JATL.EQ.IATL) GOTO 12
          XR(ND)=SYM_AT(1,JATL)-SYM_AT(1,IATL)
          YR(ND)=SYM_AT(2,JATL)-SYM_AT(2,IATL)
          ZR(ND)=SYM_AT(3,JATL)-SYM_AT(3,IATL)
          R(ND)=SQRT(XR(ND)*XR(ND)+YR(ND)*YR(ND)+ZR(ND)*ZR(ND))
          DIJ=DIJ+R(ND)
          IF((ILE.EQ.1).AND.(DIJ.GT.RLE)) IT(ND-1)=1
          IF((IT(ND-1).EQ.1).AND.(ND.GT.1)) GOTO 42
          JPOS(ND,1)=JTYP
          JPOS(ND,2)=JNUM
          JPOS(ND,3)=JATL
          NPATH(ND)=NPATH(ND)+1.
          IF(ND.GT.1) THEN
            COSTHMIJ=(XR(ND)*XR(ND-1)+YR(ND)*YR(ND-1)+ZR(ND)*ZR(ND-1))/(
     &R(ND)*R(ND-1))
            IF(IFWD.EQ.1) THEN
              IF(IBWD(ITYP).EQ.0) THEN
                IF(COSTHMIJ.LT.COSFWDI) THEN
                  NTHOF=NTHOF+1
                  IN(ND-1)=1
                  IF(NTHOF.GT.NTHOUT) THEN
                    IT(ND-1)=1
                  ENDIF
                ENDIF
              ELSEIF(IBWD(ITYP).EQ.1) THEN
                IF((COSTHMIJ.GT.COSBWDI).AND.(COSTHMIJ.LT.-SMALL)) THEN
                  NTHOF=NTHOF+1
                  IN(ND-1)=1
                  IF(NTHOF.GT.NTHOUT) THEN
                    IT(ND-1)=1
                  ENDIF
                ENDIF
                IF((COSTHMIJ.LT.COSFWDI).AND.(COSTHMIJ.GE.0.)) THEN
                  NTHOF=NTHOF+1
                  IN(ND-1)=1
                  IF(NTHOF.GT.NTHOUT) THEN
                    IT(ND-1)=1
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
          IF((IT(ND-1).EQ.1).AND.(ND.GT.1)) GOTO 42
          RHOIJ=VK(JE)*R(ND)
          CTROIS1=ZR(ND)/R(ND)
          IF(CTROIS1.GT.1.)  THEN
            CTROIS1=1.
          ELSEIF(CTROIS1.LT.-1.) THEN
            CTROIS1=-1.
          ENDIF
          THIJ=ACOS(CTROIS1)
          COMPL1= XR(ND)+IC*YR(ND)
          CALL ARCSIN(COMPL1,CTROIS1,PHIIJ)
          IF((ND.GT.1).AND.((ND-1).LT.NDIF)) THEN
            IF((IDWSPH.EQ.1).AND.(ISPEED.EQ.1)) GOTO 40
            ZSURFI=ZSURF-ZR(ND-1)
            IF(IDCM.EQ.1) THEN
              UJ2(ITYP)=SIG2(R(ND-1),ITYP)
            ENDIF
            IF(ABS(ZSURFI).LE.SMALL) THEN
              IF(ABS(COSTHMIJ-1.).GT.SMALL) THEN
                CSKZ2I=(CTROIS1-COS(THMI))*(CTROIS1-COS(THMI))/(2.-2.*CO
     &STHMIJ)
              ELSE
                CSKZ2I=1.
              ENDIF
              UII=UJ2(ITYP)*(1.+CSKZ2I*(RSJ-1.))
            ELSE
              UII=UJ2(ITYP)
            ENDIF
            IF((ISPEED.EQ.0).AND.(IDWSPH.EQ.1)) THEN
              XK2UI2=VK2(JE)*UII
              CALL DWSPH(ITYP,JE,XK2UI2,TLT,ISPEED)
            ENDIF
 40         IF(IDWSPH.EQ.1) THEN
              DW(ND-1)=1.
            ELSE
              DW(ND-1)=EXP(-VK2(JE)*UII*(1.-COSTHMIJ))
            ENDIF
          ENDIF
          IF(ND.EQ.1) THEN
            RHO01=RHOIJ
            TH01=THIJ
            PHI01=PHIIJ
            CALL DJMN2(TH01,RLM01,LF2,2)
            GOTO 30
          ENDIF
          IF(IPW.EQ.1) THEN
            CALL FACDIF(COSTHMIJ,JPOS(ND-1,1),JE,FTHETA)
            PWI=FTHETA*DW(ND-1)/R(ND)
            PW(ND)=PW(ND-1)*PWI
            CTL2=PI4*PW(ND)*CEX(1)/VK(JE)
            LMJ=LMAX(ITYP,JE)
            IF(ND.GT.NCUT) THEN
              IT(ND)=1
            ELSE
              IT(ND)=0
            ENDIF
            CALL HARSPH2(NL_M,TH01,PHI01,YLM1,LF2)
            CALL HARSPH2(NL_M,THIJ,PHIIJ,YLM2,LMJ)
            XMAXT=0.
            DO LJ=0,LMJ
              CTL=CTL2*TL(LJ,1,JTYP,JE)*YLM2(LJ,0)
              DO LF=LF1,LF2,ISTEP_LF
                PW1=CTL*YLM1(LF,0)*TL(LF,1,1,JE)
                XMAXT=AMAX1(XMAXT,CABS(PW1))
              ENDDO
            ENDDO
            IF((PCTINT*FREF-XMAXT.LT.-XCOMP).AND.(ND.GT.NCUT))IT(ND)=0
          ENDIF
          IF((IT(ND-1).EQ.1).OR.(IT(ND).EQ.1)) GOTO 42
          CALL EULER(THIJ,PHIIJ,THMI,PHIMI,AMIJ,BMIJ,CMIJ,IEULER)
          IF((I_CP.EQ.1).AND.(ND.EQ.NDIF)) I_ABS=2
          CALL MATDIF(NO,ND-1,LF2,ITYP,JTYP,JE,I_ABS,ISPEED,ISPHER,AMIJ,
     &BMIJ,CMIJ,RHOMI,RHOIJ)
  30      CEX(ND)=CEXP(IC*RHOIJ)/R(ND)
          CEXDW(ND)=CEX(ND)*DW(ND-1)
          IF((IJ.EQ.1).OR.(ND.EQ.NCUT)) THEN
            IF((I_CP.EQ.0).OR.(JATL.EQ.1)) THEN
              CALL PATHOP(JPOS,ND,JE,I_CP,RHO01,PHI01,RHOIJ,THIJ,PHIIJ,F
     &REF,IJ,DIJ,TAU)
              NPATH2(ND)=NPATH2(ND)+1.
            ENDIF
          ENDIF
          IF(ND.EQ.NDIF) GOTO 42
          I_ABS=0
C
          IF((I_CP.EQ.0).OR.(ND.NE.NDIF-1)) THEN
            N_TYP=N_PROT
          ELSE
            N_TYP=1
          ENDIF
C
          DO KTYP=1,N_TYP
            ND=ND+1
            IF(ND.GT.NDIF) GOTO 20
C
            IF((I_CP.EQ.0).OR.(ND.NE.NDIF)) THEN
              NBTYPK=NATYP(KTYP)
            ELSE
              NBTYPK=1
            ENDIF
C
            DO KNUM=1,NBTYPK
              KATL=NCORR(KNUM,KTYP)
              IF(KATL.EQ.JATL) GOTO 22
              JPOS(ND,1)=KTYP
              JPOS(ND,2)=KNUM
              JPOS(ND,3)=KATL
              XR(ND)=SYM_AT(1,KATL)-SYM_AT(1,JATL)
              YR(ND)=SYM_AT(2,KATL)-SYM_AT(2,JATL)
              ZR(ND)=SYM_AT(3,KATL)-SYM_AT(3,JATL)
              R(ND)=SQRT(XR(ND)*XR(ND)+YR(ND)*YR(ND)+ZR(ND)*ZR(ND))
              DIJ=DIJ+R(ND)
              IF((ILE.EQ.1).AND.(DIJ.GT.RLE)) IT(ND-1)=1
              IF(IT(ND-1).EQ.1) GOTO 32
              RHOJK=R(ND)*VK(JE)
              NPATH(ND)=NPATH(ND)+1.
              COSTHIJK=(XR(ND)*XR(ND-1)+YR(ND)*YR(ND-1)+ZR(ND)*ZR(ND-1))
     &/(R(ND)*R(ND-1))
              IF(IFWD.EQ.1) THEN
                IF(IBWD(JTYP).EQ.0) THEN
                  IF(COSTHIJK.LT.COSFWDJ) THEN
                    NTHOF=NTHOF+1
                    IN(ND-1)=1
                    IF(NTHOF.GT.NTHOUT) THEN
                      IT(ND-1)=1
                    ENDIF
                  ENDIF
                ELSEIF(IBWD(JTYP).EQ.1) THEN
                  IF((COSTHIJK.GT.COSBWDJ).AND. (COSTHIJK.LT.-SMALL)) TH
     &EN
                    NTHOF=NTHOF+1
                    IN(ND-1)=1
                    IF(NTHOF.GT.NTHOUT) THEN
                      IT(ND-1)=1
                    ENDIF
                  ENDIF
                  IF((COSTHIJK.LT.COSFWDJ).AND.(COSTHIJK.GE.0.))THEN
                    NTHOF=NTHOF+1
                    IN(ND-1)=1
                    IF(NTHOF.GT.NTHOUT) THEN
                      IT(ND-1)=1
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
              IF(IT(ND-1).EQ.1) GOTO 32
              CTROIS1=ZR(ND)/R(ND)
              IF(CTROIS1.GT.1)  THEN
                CTROIS1=1.
              ELSEIF(CTROIS1.LT.-1.) THEN
                CTROIS1=-1.
              ENDIF
              THJK=ACOS(CTROIS1)
              COMPL1= XR(ND)+IC*YR(ND)
              IF(ND-1.LT.NDIF) THEN
                IF((IDWSPH.EQ.1).AND.(ISPEED.EQ.1)) GOTO 50
                ZSURFJ=ZSURF-ZR(ND-1)
                IF(IDCM.EQ.1) THEN
                  UJ2(JTYP)=SIG2(R(ND-1),JTYP)
                ENDIF
                IF(ABS(ZSURFJ).LE.SMALL) THEN
                  IF(ABS(COSTHIJK-1.).GT.SMALL) THEN
                    CSKZ2J=(CTROIS1-COS(THIJ))*(CTROIS1-COS(THIJ))/(2.-2
     &.*COSTHIJK)
                  ELSE
                    CSKZ2J=1.
                  ENDIF
                  UJJ=UJ2(JTYP)*(1.+CSKZ2J*(RSJ-1.))
                ELSE
                  UJJ=UJ2(JTYP)
                ENDIF
                IF((ISPEED.EQ.0).AND.(IDWSPH.EQ.1)) THEN
                  XK2UJ2=VK2(JE)*UJJ
                  CALL DWSPH(JTYP,JE,XK2UJ2,TLT,ISPEED)
                ENDIF
 50             IF(IDWSPH.EQ.1) THEN
                  DW(ND-1)=1.
                ELSE
                  DW(ND-1)=EXP(-VK2(JE)*UJJ*(1.-COSTHIJK))
                ENDIF
              ENDIF
              IF(IPW.EQ.1) THEN
                CALL FACDIF(COSTHIJK,JPOS(ND-1,1),JE,FTHETA)
                PWI=FTHETA*DW(ND-1)/R(ND)
                PW(ND)=PW(ND-1)*PWI
                CTL2=PI4*PW(ND)*CEX(1)/VK(JE)
                CALL ARCSIN(COMPL1,CTROIS1,PHIJK)
                LMJ=LMAX(JTYP,JE)
                IF(ND.GT.NCUT) THEN
                  IT(ND)=1
                ELSE
                  IT(ND)=0
                ENDIF
                CALL HARSPH2(NL_M,TH01,PHI01,YLM1,LF2)
                CALL HARSPH2(NL_M,THJK,PHIJK,YLM2,LMJ)
                XMAXT=0.
                DO LJ=0,LMJ
                  CTL=CTL2*TL(LJ,1,KTYP,JE)*YLM2(LJ,0)
                  DO LF=LF1,LF2,ISTEP_LF
                    PW1=CTL*YLM1(LF,0)*TL(LF,1,1,JE)
                    XMAXT=AMAX1(XMAXT,CABS(PW1))
                  ENDDO
                ENDDO
                IF((PCTINT*FREF-XMAXT.LT.-XCOMP).AND.(ND.GT.NCUT))
     &              IT(ND)=0
              ENDIF
              IF((IT(ND-1).EQ.1).OR.(IT(ND).EQ.1)) GOTO 32
              IF((ND.LT.NDIF).OR.(IPW.EQ.0)) THEN
                CALL ARCSIN(COMPL1,CTROIS1,PHIJK)
              ENDIF
              CALL EULER(THJK,PHIJK,THIJ,PHIIJ,AIJK,BIJK,CIJK,IEULER)
              IF((I_CP.EQ.1).AND.(ND.EQ.NDIF)) I_ABS=2
              CALL MATDIF(NO,ND-1,LF2,JTYP,KTYP,JE,I_ABS,ISPEED,ISPHER,A
     &IJK,BIJK,CIJK,RHOIJ,RHOJK)
              CEX(ND)=CEXP(IC*RHOJK)/R(ND)
              CEXDW(ND)=CEX(ND)*DW(ND-1)
              IF((IJ.EQ.1).OR.(ND.EQ.NCUT)) THEN
                IF((I_CP.EQ.0).OR.(KATL.EQ.1)) THEN
                  CALL PATHOP(JPOS,ND,JE,I_CP,RHO01,PHI01,RHOJK,THJK,PHI
     &JK,FREF,IJ,DIJ,TAU)
                  NPATH2(ND)=NPATH2(ND)+1.
                ENDIF
              ENDIF
              IF(ND.EQ.NDIF) GOTO 32
              CALL FINDPATHS4(ND,KTYP,KATL,I_CP,R,XR,YR,ZR,RHOJK,THJK,PH
     &IJK,ZSURF,JPOS,PW,JE,FREF,DIJ,TAU)
  32          DIJ=DIJ-R(ND)
  22          IF(IN(ND-1).EQ.1) NTHOF=NTHOF-1
              IT(ND-1)=0
              IN(ND-1)=0
            ENDDO
  20        CONTINUE
            ND=ND-1
          ENDDO
  42      DIJ=DIJ-R(ND)
  12      IF(ND.GT.1) THEN
            IF(IN(ND-1).EQ.1) NTHOF=NTHOF-1
            IT(ND-1)=0
            IN(ND-1)=0
          ENDIF
       ENDDO
       ND=ND-1
      ENDDO
C
      RETURN
C
      END
