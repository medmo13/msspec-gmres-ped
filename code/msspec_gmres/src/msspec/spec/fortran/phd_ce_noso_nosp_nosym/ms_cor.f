C
C
C======================================================================
C
      SUBROUTINE MS_COR(JE,TAU)
C
C
C   This subroutine calculates the scattering path operator by 
C    the correlation expansion method.
C
C   The scattering path operator matrix of each small atom group 
C    is obtained by using LU decomposition method.
C
C   The running time of matrix inversion subroutine used in this program
C    scales with N^3, the memory occupied scales with N^2. We advise user to
C    use full MS method to get the scattering path operator, i.e. directly
C    with matrix inversion method if NGR is larger than 3. If NGR is less
C    than 4 (i.e <=3) this subroutine will gain time.
C
C   This subroutine never gain memory comparing to the subrourine INV_MAT_MS
C    as I use three large matrices stored in common, each matrix is larger or
C    as large as the matrix used in INV_MAT_MS.
C
C   As I don't find a good way to solve the group problem, where all the contribution
C    of group IGR<=NGR are collected and each small contribution has to be stored
C    for the further larger-atom-group contribution, it's better that users change the
C    parameter NGR_M which is set in included file 'spec.inc' to be NGR or NGR+1
C    where NGR is the cut-off.user insterested. this subrouitne works for NGR is less
C    than 6(<=5), if users want to calculate larger NGR, they should modify the code here
C    to make them workable, the code is marked by 'C' in each lines (about 300 lines
C    below here), users just release them until to the desired cut-off, the maximum is
C    9, however, users can enlarge it if they want to. Warning ! NGR_M set in 
C    included file should be larger than NGR and the figure listed below, don't forget
C    to compile the code after modification.
C
C   Users can modify the code to make it less memory-occupied, however, no matter they
C    do, the memories that used are more than full MS method used, so the only advantage
C    that this code has is to gain time when NGR<=3, with command 'common' used here,
C    the code will run faster.
C
C                                                      H.-F. Zhao : 2007
C
C                           (Photoelectron case)
C
C                                         Last modified : 31 Jan 2008
C
C
C
      USE DIM_MOD
      USE COOR_MOD
      USE INIT_L_MOD
      USE TRANS_MOD
      USE APPROX_MOD
      USE CORREXP_MOD
      USE Q_ARRAY_MOD
C
      COMPLEX*16 TAU1(LINMAX,LINFMAX,NATCLU_M),ONEC,ZEROC
C
      INTEGER NLM(NGR_M),ITYP(NGR_M),IGS(NGR_M)
C
      COMPLEX TAU(LINMAX,LINFMAX,NATCLU_M),TLJ
C
C
      ONEC=(1.D0,0.D0)
      ZEROC=(0.D0,0.D0)
C
      LM0=LMAX(1,JE)
      LM0=MIN(LM0,LF2)
      NRHS=(LM0+1)*(LM0+1)
C
      NGR_MAX=NGR_M
      NGR=NDIF
C
      IF(NGR_M.GT.NATCLU) THEN
        WRITE(6,*) '     --->  NGR_M should be smaller than NATCLU'
        WRITE(6,*) '     --->  it is reduced to NATCLU=',NATCLU
        NGR_MAX=NATCLU
      ENDIF
C
      IF(NGR.LT.1) THEN
        WRITE(6,*) '     --->  NGR < 1, no expansion is done'
        STOP
      ELSE
        IF(NGR.GT.NGR_MAX) THEN
          WRITE(6,*) '     --->  NGR is too large, reduce to NGR_M=',
     &                      NGR_MAX
          NGR=NGR_MAX
        ENDIF
      ENDIF
C
C      Case NGR = 1
C
      IF(NGR.EQ.1) THEN
        DO LJ=0,LM0
          ILJ=LJ*LJ+LJ+1
          TLJ=TL(LJ,1,1,JE)
          DO MJ=-LJ,LJ
            INDJ=ILJ+MJ
            TAU(INDJ,INDJ,1)=TLJ
          ENDDO
        ENDDO
C
        GOTO 100
      ENDIF
C
C     NGR >=2 case
C
C
      DO INDJ=1,NRHS
        TAU1(INDJ,INDJ,1)=DBLE(Q(1))*ONEC
      ENDDO
C
C     Constructs the group matrix and inverses it
C
      IGR=1
      LMJ=LMAX(1,JE)
      NLM(IGR)=LMJ
      INDJM=(LMJ+1)*(LMJ+1)
      ITYP(IGR)=1
      IGS(IGR)=1
C
      DO I=1,INDJM
        DO J=1,INDJM
          IF (J.EQ.I) THEN
            A(J,I)=ONEC
          ELSE
            A(J,I)=ZEROC
          ENDIF
        ENDDO
      ENDDO
C
      IGR=IGR+1
      CALL COREXP_SAVM(JE,IGR,NGR,NLM,ITYP,IGS,TAU1)
      IGR=IGR-1
C
C     TAU=TAU*tj
C
      DO KTYP=1,N_PROT
        NBTYPK=NATYP(KTYP)
        LMK=LMAX(KTYP,JE)
        INDKM=(LMK+1)*(LMK+1)
        DO KNUM=1,NBTYPK
          KATL=NCORR(KNUM,KTYP)
C
          DO LJ=0,LM0
            ILJ=LJ*LJ+LJ+1
            TLJ=TL(LJ,1,1,JE)
            DO MJ=-LJ,LJ
              INDJ=ILJ+MJ
C
              DO INDK=1,INDKM
                TAU(INDK,INDJ,KATL)=CMPLX(TAU1(INDK,INDJ,KATL))*TLJ
              ENDDO
C
            ENDDO
          ENDDO
C
        ENDDO
      ENDDO
C
 100  CONTINUE
C
      RETURN
C
      END
