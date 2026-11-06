C
C=======================================================================
C
      SUBROUTINE ACC_CONV(N_CALL,J_INI,J_FIN,N_ACC,K_ACC,K_MAX,I_CONV,
     &R_SN,SN,METH)
C
C  Use of the various acceleration scheme on a scalar series
C
C                 AITK : Aitken, modified Aitken
C                 RICH : Richardson
C                 SALZ : Salzer
C                 EPSI : epsilon
C                 EPSG : generalized epsilon
C                 RHOA : rho
C                 THET : theta
C                 LEGE : Legendre-Toeplitz
C                 CHEB : Chebyshev-Toeplitz
C                 OVER : Overholt
C                 DURB : Durbin
C                 DLEV : Levin d-transform  !!  differ by the
C                 TLEV : Levin t-transform  !!  choice of w_n
C                 ULEV : Levin u-transform  !!  the remainder
C                 VLEV : Levin v-transform  !!  estimate given
C                 ELEV : Levin e-transform  !!  by I_WN
C                 EULE : generalized Euler transform
C                 GBWT : Germain-Bronne-Wimp transform
C                 VARI : various algorithms : Weniger, BDG, iterated rho
C                 ITHE : iterated theta : Lubkin, iterated theta, modified Aitken 2
C                 EALG : E-algorithm
C
C  SN : series to be accelerated
C  XN : auxilliary series (chosen with I_XN) : interpolation points
C  GN : auxilliary series (chosen with I_GN) ---> E algorithm
C
C  N_DIV : number of SN(N+M,L), M > 0 and L < K+1,  necessary to compute S(N,K+1)
C          example: iterated Aitken method, we need SN(N,K), SN(N+1,K) and SN(N+2,K)
C                   so N_DIV=2
C          Because of this only (N_TABLE/N_DIV,N_TABLE/N_DIV) (n,k) tables are
C          meaningful. Hence the K_MAX size
C
C  COL_M : type of columns meaningful in the (n,k) table (example: even columns
C          for the epsilon algorithm)
C
C
C  Author : D. Sebilleau
C
C                                 Last modified :  1 Mar 2013
C
C
C      INCLUDE 'spec.inc'
      USE DIM_MOD
      USE CONVACC_MOD, L => LEVIN
      USE CONVTYP_MOD
C
      PARAMETER (N_METH=24)
C
      COMPLEX*16 SN(-1:N_ORD_M,-1:N_ORD_M)
C
      COMPLEX*16 ZEROC,ONEC
C
      REAL*8 R2_SN(N_ORD_M),RHO
C
      REAL R_SN(N_ORD_M)
C
      INTEGER N_DIV(N_METH)
C
      CHARACTER*3 COL_M(N_METH)
      CHARACTER*4 SCHE(N_METH)
      CHARACTER*10 NAME(N_METH),METH
C
C
      DATA SCHE /'AITK','RICH','EPSI','RHOA','THET','LEGE','CHEB',
     &           'OVER','DURB','DLEV','TLEV','ULEV','VLEV','ELEV',
     &           'EULE','GBWT','EALG','SALZ','VARI','ITHE','EPSG',
     &           'NONE','NONE','NONE'/
      DATA NAME /'  AITKEN  ','RICHARDSON','  EPSILON ','   RHO    ',
     &           '  THETA   ',' LEGENDRE ','CHEBYSHEV ',' OVERHOLT ',
     &           '  DURBIN  ',' D LEVIN  ',' T LEVIN  ',' U LEVIN  ',
     &           ' V LEVIN  ',' E LEVIN  ','  EULER   ','   GBW    ',
     &           '    E     ','  SALZER  ','  VARIA   ','ITER THETA',
     &           ' EPSILON G','          ','          ','          '/
      DATA COL_M /'ALL','ALL','EVE','EVE','EVE','ALL','ALL','ALL',
     &            'ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL',
     &            'ALL','ALL','ALL','ALL','EVE','ALL','ALL','ALL'/
      DATA N_DIV /2,1,1,2,2,1,1,1,4,1,1,1,1,1,1,1,1,1,2,4,1,1,1,1/
C
      J_NAME=0
      I_COL=0
C
      ZEROC=(0.D0,0.D0)
      ONEC=(1.D0,0.D0)
C
      DO J=J_INI,J_FIN
        R2_SN(J)=DBLE(R_SN(J))
      ENDDO
C
C  Finding the name of the method
C
      DO JM=1,N_METH
C
        IF(METHOD.EQ.SCHE(JM)) THEN
           J_NAME=JM
           K_MAX=N_MAX/N_DIV(JM)
           N_COEF=N_DIV(JM)
           IF(COL_M(JM).EQ.'EVE') I_COL=2
        ENDIF
C
      ENDDO
C
C  Initialization of array SN
C
      DO N=J_INI,J_FIN
C
        DO K=-1,J_FIN
C
         SN(N,K)=ZEROC
C
        ENDDO
C
      ENDDO
C
C  Initialisation of the schemes :
C
C   --  SN(N,0) is the series to be accelerated  --
C
      SN(0,-1)=ZEROC
      SN(0,0)=ONEC
C
      DO N=J_INI,J_FIN
C
        SN(N,-1)=ZEROC
        SN(N,0)=R2_SN(N)
C
      ENDDO
C
      CALL ACC_SCAL(J_INI,J_FIN,METHOD,SN)
C
C  Check for convergence : all results equal within ACC
C     in a N_TABLE x N_TABLE square
C
      IF(I_CONV.EQ.0) THEN
        CALL CHECK_CONV(J_INI,J_FIN,N_TABLE,I_CONV,N_ACC,K_ACC,K_MAX,
     &  N_COEF,I_COL,ACC,RHO,SN)
      ENDIF
C
      IF(METHOD(2:4).EQ.'LEV') THEN
        METH=NAME(J_NAME)(1:9)//CHAR(48+I_WN)
      ELSE
        METH=NAME(J_NAME)
      ENDIF
C
C  Incrementation of the number of calls to this subroutine
C         if convergence has not been achieved
C
      IF((I_CONV.EQ.0).OR.(ABS(I_ACC).EQ.2)) THEN
        N_CALL=N_CALL+1
      ELSE
        GOTO 10
      ENDIF
C
C  Printing the results in the check file
C
  15  FORMAT(2X,'*',3X,I3,5X,'(',E12.6,',',E12.6,')',3X,'(',E12.6,',',
     &E12.6,')',5X,'*')
  16  FORMAT(2X,'*',3X,I3,5X,'(',E12.6,',',E12.6,')',35X,'*')
  17  FORMAT(2X,'*',3X,I3,' --->','(',E12.6,',',E12.6,')',35X,'*')
  18  FORMAT(2X,'*',3X,I3,' --->','(',E12.6,',',E12.6,')',3X,'(',E12.6,
     &',',E12.6,')',5X,'*')
  19  FORMAT(2X,'*',3X,I3,25X,'(',E12.6,',',E12.6,')',20X,'*')
  21  FORMAT(2X,'*',3X,I3,25X,'(',E12.6,',',E12.6,')','  <--- 
     &convergence',2X,'*')
  25  FORMAT(2X,'*',3X,I3,5X,'(',E12.6,',',E12.6,')',3X,'(',E12.6,',',
     &E12.6,')',3X,'  *  <--- convergence')
  35  FORMAT(2X,'***************************************','************
     &************************')
  45  FORMAT(2X,'*                                      ','            
     &                       *')
  65  FORMAT(2X,'*              Exact result S = (',E12.6,',',E12.6,') 
     &              *')
  75  FORMAT(2X,'*   Order              Taylor ',20X,A10,'             
     & *')
 105  FORMAT(2X,'*    Convergence          ',A4,26X,A4,14X,'*',/2X,'*  
     &     order ',60X,'*')
 133  FORMAT(//,5X,'<<<<<<<<<<  THIS METHOD IS NOT IMPLEMENTED ',
     &'>>>>>>>>>>',//)
C
  10  RETURN
C
      END
