C
C=======================================================================
C
      SUBROUTINE HEADERS(IUO2)
C
C   This subroutine writes headers containing the main parameters
C       of the calculation in the result file. The number of
C       lines written depends of the spectroscopy
C
C                                          Last modified : 31 Jan 2013
C
      USE DIM_MOD
C
      USE ALGORITHM_MOD
      USE APPROX_MOD
      USE EXAFS_MOD
      USE HEADER_MOD
      USE INFILES_MOD
      USE INIT_A_MOD
      USE INIT_J_MOD
      USE INIT_L_MOD
      USE INIT_M_MOD
      USE MOYEN_MOD
      USE PARCAL_MOD
      USE PARCAL_A_MOD
      USE TYPCAL_MOD
      USE TYPCAL_A_MOD
      USE TYPEXP_MOD
      USE VALIN_MOD
      USE VALIN_AV_MOD
      USE VALFIN_MOD
      USE VALEX_A_MOD
C
C
C
C
C
C
C
C
C
C
C
      WRITE(IUO2,1)
      WRITE(IUO2,2)
C
C   Input files section:
C
C      Checking the size of filenames
C
      N_CHAR1=0
      DO J_CHAR=1,24
        IF(INFILE1(J_CHAR:J_CHAR).EQ.' ') GOTO 500
        N_CHAR1=N_CHAR1+1
      ENDDO
 500  CONTINUE
C
      N_CHAR2=0
      DO J_CHAR=1,24
        IF(INFILE2(J_CHAR:J_CHAR).EQ.' ') GOTO 501
        N_CHAR2=N_CHAR2+1
      ENDDO
 501  CONTINUE
C
      N_CHAR3=0
      DO J_CHAR=1,24
        IF(INFILE3(J_CHAR:J_CHAR).EQ.' ') GOTO 502
        N_CHAR3=N_CHAR3+1
      ENDDO
 502  CONTINUE
C
      N_CHAR4=0
      DO J_CHAR=1,24
        IF(INFILE4(J_CHAR:J_CHAR).EQ.' ') GOTO 503
        N_CHAR4=N_CHAR4+1
      ENDDO
 503  CONTINUE
C
      WRITE(IUO2,3) INFILE1(6:N_CHAR1)
      WRITE(IUO2,4) INFILE2(4:N_CHAR2)
      IF(INTERACT.NE.'NOINTER') THEN
        WRITE(IUO2,5) INFILE3(5:N_CHAR3)
      ENDIF
      WRITE(IUO2,6) INFILE4(6:N_CHAR4)
      WRITE(IUO2,2)
C
C   Type of calculation
C
      WRITE(IUO2,2)
C
      IF(SPECTRO.EQ.'PHD') THEN
        WRITE(IUO2,11) SPECTRO,ALGO1
        IF(ALGO1.EQ.'SE') THEN
          WRITE(IUO2,12) NO,NDIF,IFWD,IPW,ILENGTH
        ELSEIF(ALGO1.EQ.'CE') THEN
          WRITE(IUO2,13) NDIF
        ENDIF
        WRITE(IUO2,14) VINT
      ELSEIF(SPECTRO.EQ.'XAS') THEN
        WRITE(IUO2,11) SPECTRO,ALGO1
        IF(ALGO1.EQ.'SE') THEN
          WRITE(IUO2,12) NO,NDIF,IFWD,IPW,ILENGTH
        ELSEIF(ALGO1.EQ.'CE') THEN
          WRITE(IUO2,13) NDIF
        ENDIF
        WRITE(IUO2,14) VINT
      ELSEIF(SPECTRO.EQ.'LED') THEN
        WRITE(IUO2,11) SPECTRO,ALGO1
        IF(ALGO1.EQ.'SE') THEN
          WRITE(IUO2,12) NO,NDIF,IFWD,IPW,ILENGTH
        ELSEIF(ALGO1.EQ.'CE') THEN
          WRITE(IUO2,13) NDIF
        ENDIF
        WRITE(IUO2,14) VINT
      ELSEIF(SPECTRO.EQ.'AED') THEN
        WRITE(IUO2,11) SPECTRO,ALGO2
        IF(ALGO1.EQ.'SE') THEN
          WRITE(IUO2,12) NO,NDIF,IFWD,IPW,ILENGTH
        ELSEIF(ALGO1.EQ.'CE') THEN
          WRITE(IUO2,13) NDIF
        ENDIF
        WRITE(IUO2,14) VINT
      ELSEIF(SPECTRO.EQ.'APC') THEN
        WRITE(IUO2,15) SPECTRO,ALGO1,ALGO2
        WRITE(IUO2,14) VINT
      ELSEIF(SPECTRO.EQ.'EIG') THEN
        WRITE(IUO2,11) SPECTRO,ALGO1
      ELSEIF(SPECTRO.EQ.'RES') THEN
        CONTINUE
      ELSEIF(SPECTRO.EQ.'ELS') THEN
        CONTINUE
      ENDIF
C
      WRITE(IUO2,2)
C
C   Initial state parameters
C
      IF(SPECTRO.EQ.'PHD') THEN
        WRITE(IUO2,21) NI,NLI,S_O,INITL
      ELSEIF(SPECTRO.EQ.'XAS') THEN
        WRITE(IUO2,22) EDGE,NEDGE,INITL
      ELSEIF(SPECTRO.EQ.'LED') THEN
        CONTINUE
      ELSEIF(SPECTRO.EQ.'AED') THEN
        WRITE(IUO2,24) AUGER,MULTIPLET
      ELSEIF(SPECTRO.EQ.'APC') THEN
        WRITE(IUO2,21) NI,NLI,S_O,INITL
        WRITE(IUO2,24) AUGER,MULTIPLET
      ELSEIF(SPECTRO.EQ.'RES') THEN
        CONTINUE
      ELSEIF(SPECTRO.EQ.'ELS') THEN
        CONTINUE
      ENDIF
C
      WRITE(IUO2,2)
C
C   Angular and energy parameters
C
      IF(SPECTRO.EQ.'PHD') THEN
        WRITE(IUO2,35)
        WRITE(IUO2,34) THLUM,PHILUM,ELUM
        WRITE(IUO2,2)
        WRITE(IUO2,36)
        WRITE(IUO2,31) THETA0,THETA1
        WRITE(IUO2,32) PHI0,PHI1
        WRITE(IUO2,33) E0,EFIN
      ELSEIF(SPECTRO.EQ.'XAS') THEN
        WRITE(IUO2,35)
        WRITE(IUO2,33) EK_INI,EK_FIN
        WRITE(IUO2,34) THLUM,PHILUM,ELUM
      ELSEIF(SPECTRO.EQ.'LED') THEN
        WRITE(IUO2,35)
        WRITE(IUO2,31) THLUM,PHILUM
        WRITE(IUO2,2)
        WRITE(IUO2,36)
        WRITE(IUO2,31) THETA0,THETA1
        WRITE(IUO2,32) PHI0,PHI1
        WRITE(IUO2,2)
        WRITE(IUO2,33) E0,EFIN
      ELSEIF(SPECTRO.EQ.'AED') THEN
        WRITE(IUO2,36)
        WRITE(IUO2,31) THETA0_A,THETA1_A
        WRITE(IUO2,32) PHI0_A,PHI1_A
      ELSEIF(SPECTRO.EQ.'APC') THEN
        WRITE(IUO2,35)
        WRITE(IUO2,34) THLUM,PHILUM,ELUM
        WRITE(IUO2,2)
        WRITE(IUO2,37)
        WRITE(IUO2,31) THETA0,THETA1
        WRITE(IUO2,32) PHI0,PHI1
        WRITE(IUO2,33) E0,EFIN
        WRITE(IUO2,2)
        WRITE(IUO2,38)
        WRITE(IUO2,31) THETA0_A,THETA1_A
        WRITE(IUO2,32) PHI0_A,PHI1_A
      ELSEIF(SPECTRO.EQ.'EIG') THEN
        WRITE(IUO2,33) EK_INI,EK_FIN
      ELSEIF(SPECTRO.EQ.'RES') THEN
        CONTINUE
      ELSEIF(SPECTRO.EQ.'ELS') THEN
        CONTINUE
      ENDIF
C
C   End of headers
C
      WRITE(IUO2,2)
      WRITE(IUO2,1)
      WRITE(IUO2,39)
C
C   Formats
C
   1  FORMAT('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!','
     &!!!!!!!!!!!!!!!!')
   2  FORMAT('!',69X,'!')
   3  FORMAT('!',10X,'data file        :  ',A19,20X,'!')
   4  FORMAT('!',10X,'t-matrix file    :    ',A17,20X,'!')
   5  FORMAT('!',10X,'rad integral file: ',A20,20X,'!')
   6  FORMAT('!',10X,'cluster file     :  ',A19,20X,'!')
C
  11  FORMAT('!',10X,'spectroscopy     :      ',A3,8X,'algorithm : ',A2,
     &10X,'!')
  12  FORMAT('!',15X,'NO = ',I1,'  NDIF = ',I2,'  IFWD = ',I1,'  IPW = '
     &,I1,'  ILENGTH = ',I1,5X,'!')
  13  FORMAT('!',15X,'NDIF = ',I2,45X,'!')
  14  FORMAT('!',10X,'inner potential  :    ',F6.2,' eV',28X,'!')
  15  FORMAT('!',10X,'spectroscopy: ',A3,10X,'algorithm (photo): ',A2,11
     &X,'!',/,'!',37X,'algorithm (auger): ',A2, 11X,'!')
C
  21  FORMAT('!',10X,'initial state    : ',I1,A1,1X,A3,' selection rules
     &:',' INITL = ',I2,6X,'!')
  22  FORMAT('!',10X,'initial state    : ',A1,I1,2X,' selection rules:',
     &' INITL = ',I2,8X,'!')
  24  FORMAT('!',10X,'initial state    : ',A6,2X,' multiplet: ',A3,17X,'
     &!')
C
  31  FORMAT('!',10X,'THETA_INI: ',F8.2,6X,'THETA_FIN: ',F8.2,15X,'!')
  32  FORMAT('!',10X,'PHI_INI  : ',F8.2,6X,'PHI_FIN  : ',F8.2,15X,'!')
  33  FORMAT('!',10X,'E_INI    : ',F8.2,' eV',3X,'E_FIN    : ',F8.2,' eV
     &',12X,'!')
  34  FORMAT('!',10X,'THETA_LUM: ',F8.2,2X,'PHI_LUM: ',F8.2,2X,'E_LUM: '
     &,F8.2,' eV !')
  35  FORMAT('!',10X,'incoming beam    : ',40X,'!')
  36  FORMAT('!',10X,'outgoing beam    : ',40X,'!')
  37  FORMAT('!',10X,'photoelectron beam:',40X,'!')
  38  FORMAT('!',10X,'auger beam        :',40X,'!')
  39  FORMAT(71X)
C
      RETURN
C
      END
