C=======================================================================
      MODULE ADSORB_MOD
        IMPLICIT NONE
        INTEGER :: IADS
        INTEGER :: NATA
        INTEGER :: NADS1
        INTEGER :: NADS2
        INTEGER :: NADS3
        REAL, ALLOCATABLE, DIMENSION(:,:) :: ADS
        INTEGER :: NCOUCH
      CONTAINS
        SUBROUTINE ALLOC_ADSORB()
          USE DIM_MOD
          IF (ALLOCATED(ADS)) THEN
            DEALLOCATE(ADS)
          ENDIF
          ALLOCATE(ADS(3,900))
        END SUBROUTINE ALLOC_ADSORB
      END MODULE ADSORB_MOD

C=======================================================================
      MODULE APPROX_MOD
        IMPLICIT NONE
        INTEGER :: NDIF
        INTEGER :: NO
        INTEGER :: ISPHER
        INTEGER :: IFWD
        INTEGER :: NTHOUT
        REAL, ALLOCATABLE, DIMENSION(:) :: RTHFWD
        INTEGER, ALLOCATABLE, DIMENSION(:) :: IBWD
        REAL, ALLOCATABLE, DIMENSION(:) :: RTHBWD
        INTEGER :: IPW
        INTEGER :: NCUT
        REAL :: PCTINT
        INTEGER :: IPP
        INTEGER :: ISPEED
        INTEGER :: IATTS
        INTEGER :: ILENGTH
        REAL :: RLENGTH
      CONTAINS
        SUBROUTINE ALLOC_APPROX()
          USE DIM_MOD
          IF (ALLOCATED(RTHFWD)) THEN
            DEALLOCATE(RTHFWD)
          ENDIF
          ALLOCATE(RTHFWD(NATP_M))
          IF (ALLOCATED(IBWD)) THEN
            DEALLOCATE(IBWD)
          ENDIF
          ALLOCATE(IBWD(NATP_M))
          IF (ALLOCATED(RTHBWD)) THEN
            DEALLOCATE(RTHBWD)
          ENDIF
          ALLOCATE(RTHBWD(NATP_M))
        END SUBROUTINE ALLOC_APPROX
      END MODULE APPROX_MOD

C=======================================================================
      MODULE ATOMS_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:) :: VALZ
        REAL :: VALZ_MAX
        INTEGER, ALLOCATABLE, DIMENSION(:) :: NZAT
        INTEGER :: I_GR
        INTEGER :: I_INV
        CHARACTER*2, ALLOCATABLE, DIMENSION(:) :: CHEM
      CONTAINS
        SUBROUTINE ALLOC_ATOMS()
          USE DIM_MOD
          IF (ALLOCATED(VALZ)) THEN
            DEALLOCATE(VALZ)
          ENDIF
          ALLOCATE(VALZ(NATCLU_M))
          IF (ALLOCATED(NZAT)) THEN
            DEALLOCATE(NZAT)
          ENDIF
          ALLOCATE(NZAT(NATCLU_M))
          IF (ALLOCATED(CHEM)) THEN
            DEALLOCATE(CHEM)
          ENDIF
          ALLOCATE(CHEM(NATCLU_M))
        END SUBROUTINE ALLOC_ATOMS
      END MODULE ATOMS_MOD

C=======================================================================
      MODULE AUGER_MOD
        IMPLICIT NONE
        INTEGER, ALLOCATABLE, DIMENSION(:) :: NLIN_A
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: L_BOUNDS
        CHARACTER*6 :: AUGER
      CONTAINS
        SUBROUTINE ALLOC_AUGER()
          USE DIM_MOD
          IF (ALLOCATED(NLIN_A)) THEN
            DEALLOCATE(NLIN_A)
          ENDIF
          ALLOCATE(NLIN_A(0:20))
          IF (ALLOCATED(L_BOUNDS)) THEN
            DEALLOCATE(L_BOUNDS)
          ENDIF
          ALLOCATE(L_BOUNDS(0:20,2))
        END SUBROUTINE ALLOC_AUGER
      END MODULE AUGER_MOD

C=======================================================================
      MODULE BASES_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:) :: ATBAS
        REAL, ALLOCATABLE, DIMENSION(:) :: VECBAS
      CONTAINS
        SUBROUTINE ALLOC_BASES()
          USE DIM_MOD
          IF (ALLOCATED(ATBAS)) THEN
            DEALLOCATE(ATBAS)
          ENDIF
          ALLOCATE(ATBAS(3*NATP_M))
          IF (ALLOCATED(VECBAS)) THEN
            DEALLOCATE(VECBAS)
          ENDIF
          ALLOCATE(VECBAS(9))
        END SUBROUTINE ALLOC_BASES
      END MODULE BASES_MOD

C=======================================================================
      MODULE CLUSLIM_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:) :: X_MAX
        REAL, ALLOCATABLE, DIMENSION(:) :: X_MIN
        REAL, ALLOCATABLE, DIMENSION(:) :: Y_MAX
        REAL, ALLOCATABLE, DIMENSION(:) :: Y_MIN
        REAL, ALLOCATABLE, DIMENSION(:) :: VAL
        INTEGER :: NPLAN
      CONTAINS
        SUBROUTINE ALLOC_CLUSLIM()
          USE DIM_MOD
          IF (ALLOCATED(X_MAX)) THEN
            DEALLOCATE(X_MAX)
          ENDIF
          ALLOCATE(X_MAX(NATCLU_M))
          IF (ALLOCATED(X_MIN)) THEN
            DEALLOCATE(X_MIN)
          ENDIF
          ALLOCATE(X_MIN(NATCLU_M))
          IF (ALLOCATED(Y_MAX)) THEN
            DEALLOCATE(Y_MAX)
          ENDIF
          ALLOCATE(Y_MAX(NATCLU_M))
          IF (ALLOCATED(Y_MIN)) THEN
            DEALLOCATE(Y_MIN)
          ENDIF
          ALLOCATE(Y_MIN(NATCLU_M))
          IF (ALLOCATED(VAL)) THEN
            DEALLOCATE(VAL)
          ENDIF
          ALLOCATE(VAL(NATCLU_M))
        END SUBROUTINE ALLOC_CLUSLIM
      END MODULE CLUSLIM_MOD

C=======================================================================
      MODULE COOR_MOD
        IMPLICIT NONE
        INTEGER :: NATCLU
        INTEGER :: N_PROT
        INTEGER, ALLOCATABLE, DIMENSION(:) :: NATYP
        INTEGER, ALLOCATABLE, DIMENSION(:) :: NCHTYP
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NCORR
        INTEGER, ALLOCATABLE, DIMENSION(:) :: INEW_AT
        REAL, ALLOCATABLE, DIMENSION(:,:) :: SYM_AT
      CONTAINS
        SUBROUTINE ALLOC_COOR()
          USE DIM_MOD
          IF (ALLOCATED(NATYP)) THEN
            DEALLOCATE(NATYP)
          ENDIF
          ALLOCATE(NATYP(NATM))
          IF (ALLOCATED(NCHTYP)) THEN
            DEALLOCATE(NCHTYP)
          ENDIF
          ALLOCATE(NCHTYP(NATP_M))
          IF (ALLOCATED(NCORR)) THEN
            DEALLOCATE(NCORR)
          ENDIF
          ALLOCATE(NCORR(NAT_EQ_M,NATP_M))
          IF (ALLOCATED(INEW_AT)) THEN
            DEALLOCATE(INEW_AT)
          ENDIF
          ALLOCATE(INEW_AT(NATCLU_M))
          IF (ALLOCATED(SYM_AT)) THEN
            DEALLOCATE(SYM_AT)
          ENDIF
          ALLOCATE(SYM_AT(3,NATCLU_M))
        END SUBROUTINE ALLOC_COOR
      END MODULE COOR_MOD

C=======================================================================
      MODULE CORREXP_MOD
        IMPLICIT NONE
        COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: A
      CONTAINS
        SUBROUTINE ALLOC_CORREXP()
          USE DIM_MOD
          IF (ALLOCATED(A)) THEN
            DEALLOCATE(A)
          ENDIF
          ALLOCATE(A(NLMM,NLMM))
        END SUBROUTINE ALLOC_CORREXP
      END MODULE CORREXP_MOD

C=======================================================================
      MODULE DEBWAL_MOD
        IMPLICIT NONE
        INTEGER :: IDCM
        INTEGER :: IDWSPH
        REAL :: TD
        REAL :: QD
        REAL :: TEMP
        REAL :: RSJ
        REAL, ALLOCATABLE, DIMENSION(:) :: UJ2
      CONTAINS
        SUBROUTINE ALLOC_DEBWAL()
          USE DIM_MOD
          IF (ALLOCATED(UJ2)) THEN
            DEALLOCATE(UJ2)
          ENDIF
          ALLOCATE(UJ2(NATM))
        END SUBROUTINE ALLOC_DEBWAL
      END MODULE DEBWAL_MOD

C=======================================================================
      MODULE INDAT_MOD
        IMPLICIT NONE
        CHARACTER*24, ALLOCATABLE, DIMENSION(:) :: INDATA
      CONTAINS
        SUBROUTINE ALLOC_INDAT()
          USE DIM_MOD
          IF (ALLOCATED(INDATA)) THEN
            DEALLOCATE(INDATA)
          ENDIF
          ALLOCATE(INDATA(100))
        END SUBROUTINE ALLOC_INDAT
      END MODULE INDAT_MOD

C=======================================================================
      MODULE INIT_A_MOD
        IMPLICIT NONE
        INTEGER :: LI_C
        INTEGER :: LI_I
        INTEGER :: LI_A
      CONTAINS
        SUBROUTINE ALLOC_INIT_A()
          USE DIM_MOD
        END SUBROUTINE ALLOC_INIT_A
      END MODULE INIT_A_MOD

C=======================================================================
      MODULE INIT_L_MOD
        IMPLICIT NONE
        INTEGER :: LI
        INTEGER :: INITL
        INTEGER :: NNL
        INTEGER :: LF1
        INTEGER :: LF2
        INTEGER :: ISTEP_LF
      CONTAINS
        SUBROUTINE ALLOC_INIT_L()
          USE DIM_MOD
        END SUBROUTINE ALLOC_INIT_L
      END MODULE INIT_L_MOD

C=======================================================================
      MODULE INIT_J_MOD
        IMPLICIT NONE
        INTEGER :: JF1
        INTEGER :: JF2
        INTEGER :: I_SO
        CHARACTER*3 :: S_O
      CONTAINS
        SUBROUTINE ALLOC_INIT_J()
          USE DIM_MOD
        END SUBROUTINE ALLOC_INIT_J
      END MODULE INIT_J_MOD

C=======================================================================
      MODULE INIT_M_MOD
        IMPLICIT NONE
        INTEGER :: I_SHELL
        INTEGER :: I_MULT
        INTEGER :: L_MUL
        INTEGER :: J_MUL
        INTEGER :: S_MUL
        CHARACTER*3 :: MULTIPLET
      CONTAINS
        SUBROUTINE ALLOC_INIT_M()
          USE DIM_MOD
        END SUBROUTINE ALLOC_INIT_M
      END MODULE INIT_M_MOD

C=======================================================================
      MODULE INFILES_MOD
        IMPLICIT NONE
        CHARACTER*24 :: INFILE1
        CHARACTER*24 :: INFILE2
        CHARACTER*24 :: INFILE3
        CHARACTER*24 :: INFILE4
        CHARACTER*24 :: INFILE5
        CHARACTER*24 :: INFILE6
        CHARACTER*24 :: INFILE7
        CHARACTER*24 :: INFILE8
        CHARACTER*24 :: INFILE9
      CONTAINS
        SUBROUTINE ALLOC_INFILES()
          USE DIM_MOD
        END SUBROUTINE ALLOC_INFILES
      END MODULE INFILES_MOD

C=======================================================================
      MODULE INUNITS_MOD
        IMPLICIT NONE
        INTEGER :: IUI1
        INTEGER :: IUI2
        INTEGER :: IUI3
        INTEGER :: IUI4
        INTEGER :: IUI5
        INTEGER :: IUI6
        INTEGER :: IUI7
        INTEGER :: IUI8
        INTEGER :: IUI9
      CONTAINS
        SUBROUTINE ALLOC_INUNITS()
          USE DIM_MOD
        END SUBROUTINE ALLOC_INUNITS
      END MODULE INUNITS_MOD

C=======================================================================
      MODULE LIMAMA_MOD
        IMPLICIT NONE
        INTEGER :: NIV
        REAL :: COUPUR
      CONTAINS
        SUBROUTINE ALLOC_LIMAMA()
          USE DIM_MOD
        END SUBROUTINE ALLOC_LIMAMA
      END MODULE LIMAMA_MOD

C=======================================================================
      MODULE LPMOY_MOD
        IMPLICIT NONE
        INTEGER :: ILPM
        INTEGER :: NZA
        REAL :: XMTA
        REAL :: RHOTA
        REAL :: XLPM0
      CONTAINS
        SUBROUTINE ALLOC_LPMOY()
          USE DIM_MOD
        END SUBROUTINE ALLOC_LPMOY
      END MODULE LPMOY_MOD

C=======================================================================
      MODULE MASSAT_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:) :: XMT
      CONTAINS
        SUBROUTINE ALLOC_MASSAT()
          USE DIM_MOD
          IF (ALLOCATED(XMT)) THEN
            DEALLOCATE(XMT)
          ENDIF
          ALLOCATE(XMT(NATM))
        END SUBROUTINE ALLOC_MASSAT
      END MODULE MASSAT_MOD

C=======================================================================
      MODULE MILLER_MOD
        IMPLICIT NONE
        INTEGER :: IH
        INTEGER :: IK
        INTEGER :: II
        INTEGER :: IL
        INTEGER :: IVG0
        INTEGER, ALLOCATABLE, DIMENSION(:) :: IVN
      CONTAINS
        SUBROUTINE ALLOC_MILLER()
          USE DIM_MOD
          IF (ALLOCATED(IVN)) THEN
            DEALLOCATE(IVN)
          ENDIF
          ALLOCATE(IVN(3))
        END SUBROUTINE ALLOC_MILLER
      END MODULE MILLER_MOD

C=======================================================================
      MODULE OUTUNITS_MOD
        IMPLICIT NONE
        INTEGER :: IUO1
        INTEGER :: IUO2
        INTEGER :: IUO3
        INTEGER :: IUO4
        INTEGER :: IUSCR
        INTEGER :: IUSCR2
      CONTAINS
        SUBROUTINE ALLOC_OUTUNITS()
          USE DIM_MOD
        END SUBROUTINE ALLOC_OUTUNITS
      END MODULE OUTUNITS_MOD

C=======================================================================
      MODULE PARCAL_MOD
        IMPLICIT NONE
        INTEGER :: NPHI
        INTEGER :: NE
        INTEGER :: NTHETA
        INTEGER :: NFTHET
        INTEGER :: NEPS
      CONTAINS
        SUBROUTINE ALLOC_PARCAL()
          USE DIM_MOD
        END SUBROUTINE ALLOC_PARCAL
      END MODULE PARCAL_MOD

C=======================================================================
      MODULE PARCAL_A_MOD
        IMPLICIT NONE
        INTEGER :: NPHI_A
        INTEGER :: NE_A
        INTEGER :: NTHETA_A
        INTEGER :: NFTHET_A
      CONTAINS
        SUBROUTINE ALLOC_PARCAL_A()
          USE DIM_MOD
        END SUBROUTINE ALLOC_PARCAL_A
      END MODULE PARCAL_A_MOD

C=======================================================================
      MODULE Q_ARRAY_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:) :: Q
      CONTAINS
        SUBROUTINE ALLOC_Q_ARRAY()
          USE DIM_MOD
          IF (ALLOCATED(Q)) THEN
            DEALLOCATE(Q)
          ENDIF
          ALLOCATE(Q(NGR_M))
        END SUBROUTINE ALLOC_Q_ARRAY
      END MODULE Q_ARRAY_MOD

C=======================================================================
      MODULE RELADS_MOD
        IMPLICIT NONE
        INTEGER :: NRELA
        REAL, ALLOCATABLE, DIMENSION(:) :: PCRELA
      CONTAINS
        SUBROUTINE ALLOC_RELADS()
          USE DIM_MOD
          IF (ALLOCATED(PCRELA)) THEN
            DEALLOCATE(PCRELA)
          ENDIF
          ALLOCATE(PCRELA(3))
        END SUBROUTINE ALLOC_RELADS
      END MODULE RELADS_MOD

C=======================================================================
      MODULE RELAX_MOD
        IMPLICIT NONE
        INTEGER :: IREL
        INTEGER :: NREL
        REAL, ALLOCATABLE, DIMENSION(:) :: PCREL
        REAL :: OMEGA1
        REAL :: OMEGA2
      CONTAINS
        SUBROUTINE ALLOC_RELAX()
          USE DIM_MOD
          IF (ALLOCATED(PCREL)) THEN
            DEALLOCATE(PCREL)
          ENDIF
          ALLOCATE(PCREL(10))
        END SUBROUTINE ALLOC_RELAX
      END MODULE RELAX_MOD

C=======================================================================
      MODULE RESEAU_MOD
        IMPLICIT NONE
        INTEGER :: NCRIST
        INTEGER :: NCENTR
        INTEGER :: IBAS
        INTEGER :: NAT
        REAL :: A
        REAL :: BSURA
        REAL :: CSURA
        CHARACTER*3 :: UNIT
      CONTAINS
        SUBROUTINE ALLOC_RESEAU()
          USE DIM_MOD
        END SUBROUTINE ALLOC_RESEAU
      END MODULE RESEAU_MOD

C=======================================================================
      MODULE SPIN_MOD
        IMPLICIT NONE
        INTEGER :: ISPIN
        INTEGER :: IDICHR
        INTEGER :: NSPIN
        INTEGER :: NSPIN2
        INTEGER :: ISFLIP
        INTEGER :: IR_DIA
        INTEGER :: NSTEP
      CONTAINS
        SUBROUTINE ALLOC_SPIN()
          USE DIM_MOD
        END SUBROUTINE ALLOC_SPIN
      END MODULE SPIN_MOD

C=======================================================================
      MODULE TESTS_MOD
        IMPLICIT NONE
        INTEGER :: ITEST
        INTEGER :: IPRINT
        INTEGER :: ISORT1
        INTEGER :: NPATHP
        INTEGER :: ISOM
      CONTAINS
        SUBROUTINE ALLOC_TESTS()
          USE DIM_MOD
        END SUBROUTINE ALLOC_TESTS
      END MODULE TESTS_MOD

C=======================================================================
      MODULE TRANS_MOD
        IMPLICIT NONE
        COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: DLT
        COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: TL
        COMPLEX, ALLOCATABLE, DIMENSION(:) :: VK
        REAL, ALLOCATABLE, DIMENSION(:) :: VK2
        INTEGER :: IPOTC
        INTEGER :: ITL
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LMAX
      CONTAINS
        SUBROUTINE ALLOC_TRANS()
          USE DIM_MOD
          IF (ALLOCATED(DLT)) THEN
            DEALLOCATE(DLT)
          ENDIF
          ALLOCATE(DLT(NE_M,NATM,0:18,2))
          IF (ALLOCATED(TL)) THEN
            DEALLOCATE(TL)
          ENDIF
          ALLOCATE(TL(0:NT_M,4,NATM,NE_M))
          IF (ALLOCATED(VK)) THEN
            DEALLOCATE(VK)
          ENDIF
          ALLOCATE(VK(NE_M))
          IF (ALLOCATED(VK2)) THEN
            DEALLOCATE(VK2)
          ENDIF
          ALLOCATE(VK2(NE_M))
          IF (ALLOCATED(LMAX)) THEN
            DEALLOCATE(LMAX)
          ENDIF
          ALLOCATE(LMAX(NATM,NE_M))
        END SUBROUTINE ALLOC_TRANS
      END MODULE TRANS_MOD

C=======================================================================
      MODULE TL_AED_MOD
        IMPLICIT NONE
        COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: DLT_A
        COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: TL_A
        COMPLEX, ALLOCATABLE, DIMENSION(:) :: VK_A
        REAL, ALLOCATABLE, DIMENSION(:) :: VK2_A
        INTEGER :: IPOTC_A
        INTEGER :: ITL_A
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LMAX_A
      CONTAINS
        SUBROUTINE ALLOC_TL_AED()
          USE DIM_MOD
          IF (ALLOCATED(DLT_A)) THEN
            DEALLOCATE(DLT_A)
          ENDIF
          ALLOCATE(DLT_A(NE_M,NATM,0:18,2))
          IF (ALLOCATED(TL_A)) THEN
            DEALLOCATE(TL_A)
          ENDIF
          ALLOCATE(TL_A(0:NT_M,4,NATM,NE_M))
          IF (ALLOCATED(VK_A)) THEN
            DEALLOCATE(VK_A)
          ENDIF
          ALLOCATE(VK_A(NE_M))
          IF (ALLOCATED(VK2_A)) THEN
            DEALLOCATE(VK2_A)
          ENDIF
          ALLOCATE(VK2_A(NE_M))
          IF (ALLOCATED(LMAX_A)) THEN
            DEALLOCATE(LMAX_A)
          ENDIF
          ALLOCATE(LMAX_A(NATM,NE_M))
        END SUBROUTINE ALLOC_TL_AED
      END MODULE TL_AED_MOD

C=======================================================================
      MODULE TYPCAL_MOD
        IMPLICIT NONE
        INTEGER :: IPHI
        INTEGER :: IE
        INTEGER :: ITHETA
        INTEGER :: IFTHET
        INTEGER :: IMOD
        INTEGER :: IPOL
        INTEGER :: I_CP
        INTEGER :: I_EXT
        INTEGER :: I_TEST
      CONTAINS
        SUBROUTINE ALLOC_TYPCAL()
          USE DIM_MOD
        END SUBROUTINE ALLOC_TYPCAL
      END MODULE TYPCAL_MOD

C=======================================================================
      MODULE TYPCAL_A_MOD
        IMPLICIT NONE
        INTEGER :: IPHI_A
        INTEGER :: IE_A
        INTEGER :: ITHETA_A
        INTEGER :: IFTHET_A
        INTEGER :: IMOD_A
        INTEGER :: I_CP_A
        INTEGER :: I_EXT_A
        INTEGER :: I_TEST_A
      CONTAINS
        SUBROUTINE ALLOC_TYPCAL_A()
          USE DIM_MOD
        END SUBROUTINE ALLOC_TYPCAL_A
      END MODULE TYPCAL_A_MOD

C=======================================================================
      MODULE TYPEM_MOD
        IMPLICIT NONE
        INTEGER :: NEMET
        INTEGER :: IESURF
        INTEGER, ALLOCATABLE, DIMENSION(:) :: IEMET
      CONTAINS
        SUBROUTINE ALLOC_TYPEM()
          USE DIM_MOD
          IF (ALLOCATED(IEMET)) THEN
            DEALLOCATE(IEMET)
          ENDIF
          ALLOCATE(IEMET(NEMET_M))
        END SUBROUTINE ALLOC_TYPEM
      END MODULE TYPEM_MOD

C=======================================================================
      MODULE TYPEXP_MOD
        IMPLICIT NONE
        CHARACTER*3 :: SPECTRO
        CHARACTER*7 :: INTERACT
        CHARACTER*3 :: STEREO
      CONTAINS
        SUBROUTINE ALLOC_TYPEXP()
          USE DIM_MOD
        END SUBROUTINE ALLOC_TYPEXP
      END MODULE TYPEXP_MOD

C=======================================================================
      MODULE VALIN_MOD
        IMPLICIT NONE
        REAL :: PHI0
        REAL :: E0
        REAL :: THETA0
        REAL :: THLUM
        REAL :: PHILUM
        REAL :: ELUM
        REAL :: VINT
        INTEGER, ALLOCATABLE, DIMENSION(:) :: NONVOL
      CONTAINS
        SUBROUTINE ALLOC_VALIN()
          USE DIM_MOD
          IF (ALLOCATED(NONVOL)) THEN
            DEALLOCATE(NONVOL)
          ENDIF
          ALLOCATE(NONVOL(100))
        END SUBROUTINE ALLOC_VALIN
      END MODULE VALIN_MOD

C=======================================================================
      MODULE XMRHO_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:) :: XMAT
        REAL, ALLOCATABLE, DIMENSION(:) :: RHOAT
      CONTAINS
        SUBROUTINE ALLOC_XMRHO()
          USE DIM_MOD
          IF (ALLOCATED(XMAT)) THEN
            DEALLOCATE(XMAT)
          ENDIF
          ALLOCATE(XMAT(0:99))
          IF (ALLOCATED(RHOAT)) THEN
            DEALLOCATE(RHOAT)
          ENDIF
          ALLOCATE(RHOAT(0:99))
        END SUBROUTINE ALLOC_XMRHO
      END MODULE XMRHO_MOD

C=======================================================================
      MODULE CRANGL_MOD
        IMPLICIT NONE
        REAL :: ALPHAD
        REAL :: BETAD
        REAL :: GAMMAD
      CONTAINS
        SUBROUTINE ALLOC_CRANGL()
          USE DIM_MOD
        END SUBROUTINE ALLOC_CRANGL
      END MODULE CRANGL_MOD

C=======================================================================
      MODULE VECSYS_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:) :: ASYS
        REAL, ALLOCATABLE, DIMENSION(:) :: BSYS
        REAL, ALLOCATABLE, DIMENSION(:) :: CSYS
      CONTAINS
        SUBROUTINE ALLOC_VECSYS()
          USE DIM_MOD
          IF (ALLOCATED(ASYS)) THEN
            DEALLOCATE(ASYS)
          ENDIF
          ALLOCATE(ASYS(3))
          IF (ALLOCATED(BSYS)) THEN
            DEALLOCATE(BSYS)
          ENDIF
          ALLOCATE(BSYS(3))
          IF (ALLOCATED(CSYS)) THEN
            DEALLOCATE(CSYS)
          ENDIF
          ALLOCATE(CSYS(3))
        END SUBROUTINE ALLOC_VECSYS
      END MODULE VECSYS_MOD

C=======================================================================
      MODULE VIBRAT_MOD
        IMPLICIT NONE
        INTEGER, ALLOCATABLE, DIMENSION(:) :: I_FREE
      CONTAINS
        SUBROUTINE ALLOC_VIBRAT()
          USE DIM_MOD
          IF (ALLOCATED(I_FREE)) THEN
            DEALLOCATE(I_FREE)
          ENDIF
          ALLOCATE(I_FREE(NATP_M))
        END SUBROUTINE ALLOC_VIBRAT
      END MODULE VIBRAT_MOD

C=======================================================================
      MODULE COEFRLM_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: CF
      CONTAINS
        SUBROUTINE ALLOC_COEFRLM()
          USE DIM_MOD
          IF (ALLOCATED(CF)) THEN
            DEALLOCATE(CF)
          ENDIF
          ALLOCATE(CF(0:2*NL_M-2,0:2*NL_M-2,0:2*NL_M-2))
        END SUBROUTINE ALLOC_COEFRLM
      END MODULE COEFRLM_MOD

C=======================================================================
      MODULE EXPROT_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:,:) :: EXPR
      CONTAINS
        SUBROUTINE ALLOC_EXPROT()
          USE DIM_MOD
          IF (ALLOCATED(EXPR)) THEN
            DEALLOCATE(EXPR)
          ENDIF
          ALLOCATE(EXPR(0:2*NL_M-2,0:2*NL_M-2))
        END SUBROUTINE ALLOC_EXPROT
      END MODULE EXPROT_MOD

C=======================================================================
      MODULE EXPFAC_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:,:) :: EXPF
      CONTAINS
        SUBROUTINE ALLOC_EXPFAC()
          USE DIM_MOD
          IF (ALLOCATED(EXPF)) THEN
            DEALLOCATE(EXPF)
          ENDIF
          ALLOCATE(EXPF(0:2*NL_M-2,0:2*NL_M-2))
        END SUBROUTINE ALLOC_EXPFAC
      END MODULE EXPFAC_MOD

C=======================================================================
      MODULE DEXPFAC_MOD
        IMPLICIT NONE
        REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DEXPF
      CONTAINS
        SUBROUTINE ALLOC_DEXPFAC()
          USE DIM_MOD
          IF (ALLOCATED(DEXPF)) THEN
            DEALLOCATE(DEXPF)
          ENDIF
          ALLOCATE(DEXPF(0:2*NL_M-2,0:2*NL_M-2))
        END SUBROUTINE ALLOC_DEXPFAC
      END MODULE DEXPFAC_MOD

C=======================================================================
      MODULE GAUNT_C_MOD
        IMPLICIT NONE
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GNT
      CONTAINS
        SUBROUTINE ALLOC_GAUNT_C()
          USE DIM_MOD
          IF (ALLOCATED(GNT)) THEN
            DEALLOCATE(GNT)
          ENDIF
          ALLOCATE(GNT(0:N_GAUNT,LINMAX,LINMAX))
        END SUBROUTINE ALLOC_GAUNT_C
      END MODULE GAUNT_C_MOD

C=======================================================================
      MODULE LOGAMAD_MOD
        IMPLICIT NONE
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: GLD
      CONTAINS
        SUBROUTINE ALLOC_LOGAMAD()
          USE DIM_MOD
          IF (ALLOCATED(GLD)) THEN
            DEALLOCATE(GLD)
          ENDIF
          ALLOCATE(GLD(0:N_GAUNT,2))
        END SUBROUTINE ALLOC_LOGAMAD
      END MODULE LOGAMAD_MOD

C=======================================================================
      MODULE EXPFAC2_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:,:) :: EXPF2
      CONTAINS
        SUBROUTINE ALLOC_EXPFAC2()
          USE DIM_MOD
          IF (ALLOCATED(EXPF2)) THEN
            DEALLOCATE(EXPF2)
          ENDIF
          ALLOCATE(EXPF2(0:2*NL_M-2,0:2*NL_M-2))
        END SUBROUTINE ALLOC_EXPFAC2
      END MODULE EXPFAC2_MOD

C=======================================================================
      MODULE FACTSQ_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:) :: FSQ
      CONTAINS
        SUBROUTINE ALLOC_FACTSQ()
          USE DIM_MOD
          IF (ALLOCATED(FSQ)) THEN
            DEALLOCATE(FSQ)
          ENDIF
          ALLOCATE(FSQ(0:2*NL_M-2))
        END SUBROUTINE ALLOC_FACTSQ
      END MODULE FACTSQ_MOD

C=======================================================================
      MODULE ALGORITHM_MOD
        IMPLICIT NONE
        CHARACTER*2 :: ALGO1
        CHARACTER*2 :: ALGO2
        CHARACTER*2 :: ALGO3
        CHARACTER*2 :: ALGO4
      CONTAINS
        SUBROUTINE ALLOC_ALGORITHM()
          USE DIM_MOD
        END SUBROUTINE ALLOC_ALGORITHM
      END MODULE ALGORITHM_MOD

C=======================================================================
      MODULE EXAFS_MOD
        IMPLICIT NONE
        INTEGER :: NE_X
        REAL :: EK_INI
        REAL :: EK_FIN
        REAL :: EPH_INI
      CONTAINS
        SUBROUTINE ALLOC_EXAFS()
          USE DIM_MOD
        END SUBROUTINE ALLOC_EXAFS
      END MODULE EXAFS_MOD

C=======================================================================
      MODULE HEADER_MOD
        IMPLICIT NONE
        INTEGER :: NI
        CHARACTER*1 :: NLI
        CHARACTER*6 :: AUGER
        CHARACTER*1 :: EDGE
        INTEGER :: NEDGE
      CONTAINS
        SUBROUTINE ALLOC_HEADER()
          USE DIM_MOD
        END SUBROUTINE ALLOC_HEADER
      END MODULE HEADER_MOD

C=======================================================================
      MODULE MOYEN_MOD
        IMPLICIT NONE
        INTEGER :: IMOY
        INTEGER :: NDIR
        REAL :: ACCEPT
        INTEGER :: ICHKDIR
      CONTAINS
        SUBROUTINE ALLOC_MOYEN()
          USE DIM_MOD
        END SUBROUTINE ALLOC_MOYEN
      END MODULE MOYEN_MOD

C=======================================================================
      MODULE VALIN_AV_MOD
        IMPLICIT NONE
        INTEGER :: I_SET
        REAL, ALLOCATABLE, DIMENSION(:) :: TH_0
        REAL, ALLOCATABLE, DIMENSION(:) :: PH_0
      CONTAINS
        SUBROUTINE ALLOC_VALIN_AV()
          USE DIM_MOD
          IF (ALLOCATED(TH_0)) THEN
            DEALLOCATE(TH_0)
          ENDIF
          ALLOCATE(TH_0(NTH_M))
          IF (ALLOCATED(PH_0)) THEN
            DEALLOCATE(PH_0)
          ENDIF
          ALLOCATE(PH_0(NPH_M))
        END SUBROUTINE ALLOC_VALIN_AV
      END MODULE VALIN_AV_MOD

C=======================================================================
      MODULE VALFIN_MOD
        IMPLICIT NONE
        REAL :: PHI1
        REAL :: EFIN
        REAL :: THETA1
      CONTAINS
        SUBROUTINE ALLOC_VALFIN()
          USE DIM_MOD
        END SUBROUTINE ALLOC_VALFIN
      END MODULE VALFIN_MOD

C=======================================================================
      MODULE VALEX_A_MOD
        IMPLICIT NONE
        REAL :: PHI0_A
        REAL :: THETA0_A
        REAL :: PHI1_A
        REAL :: THETA1_A
      CONTAINS
        SUBROUTINE ALLOC_VALEX_A()
          USE DIM_MOD
        END SUBROUTINE ALLOC_VALEX_A
      END MODULE VALEX_A_MOD

C=======================================================================
      MODULE AMPLI_MOD
        IMPLICIT NONE
        INTEGER :: I_AMP
      CONTAINS
        SUBROUTINE ALLOC_AMPLI()
          USE DIM_MOD
        END SUBROUTINE ALLOC_AMPLI
      END MODULE AMPLI_MOD

C=======================================================================
      MODULE CONVACC_MOD
        IMPLICIT NONE
        COMPLEX*16 :: ALPHA
        COMPLEX*16 :: BETA
        INTEGER :: I_XN
        INTEGER :: I_VA
        INTEGER :: I_GN
        INTEGER :: I_WN
        INTEGER :: LEVIN
      CONTAINS
        SUBROUTINE ALLOC_CONVACC()
          USE DIM_MOD
        END SUBROUTINE ALLOC_CONVACC
      END MODULE CONVACC_MOD

C=======================================================================
      MODULE CONVTYP_MOD
        IMPLICIT NONE
        REAL*8 :: SHIFT
        REAL*8 :: ACC
        REAL*8 :: EXPO
        INTEGER :: I_PWM
        INTEGER :: I_ACC
        INTEGER :: N_ONE
        INTEGER :: N_MAX
        INTEGER :: N_ITER
        INTEGER :: N_TABLE
        CHARACTER*4 :: METHOD
      CONTAINS
        SUBROUTINE ALLOC_CONVTYP()
          USE DIM_MOD
        END SUBROUTINE ALLOC_CONVTYP
      END MODULE CONVTYP_MOD

C=======================================================================
      MODULE C_G_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: CG
      CONTAINS
        SUBROUTINE ALLOC_C_G()
          USE DIM_MOD
          IF (ALLOCATED(CG)) THEN
            DEALLOCATE(CG)
          ENDIF
          ALLOCATE(CG(-LI_M:LI_M+1,2,2))
        END SUBROUTINE ALLOC_C_G
      END MODULE C_G_MOD

C=======================================================================
      MODULE C_G_A_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: CGA
      CONTAINS
        SUBROUTINE ALLOC_C_G_A()
          USE DIM_MOD
          IF (ALLOCATED(CGA)) THEN
            DEALLOCATE(CGA)
          ENDIF
          ALLOCATE(CGA(0:NCG_M,-NCG_M:NCG_M,0:NCG_M,-NCG_M:NCG_M,0:
     &    2*NCG_M))
        END SUBROUTINE ALLOC_C_G_A
      END MODULE C_G_A_MOD

C=======================================================================
      MODULE C_G_M_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: CG_S
      CONTAINS
        SUBROUTINE ALLOC_C_G_M()
          USE DIM_MOD
          IF (ALLOCATED(CG_S)) THEN
            DEALLOCATE(CG_S)
          ENDIF
          ALLOCATE(CG_S(2,2,2))
        END SUBROUTINE ALLOC_C_G_M
      END MODULE C_G_M_MOD

C=======================================================================
      MODULE DEXPFAC2_MOD
        IMPLICIT NONE
        REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DEXPF2
      CONTAINS
        SUBROUTINE ALLOC_DEXPFAC2()
          USE DIM_MOD
          IF (ALLOCATED(DEXPF2)) THEN
            DEALLOCATE(DEXPF2)
          ENDIF
          ALLOCATE(DEXPF2(0:2*NL_M-2,0:2*NL_M-2))
        END SUBROUTINE ALLOC_DEXPFAC2
      END MODULE DEXPFAC2_MOD

C=======================================================================
      MODULE DFACTSQ_MOD
        IMPLICIT NONE
        REAL*8, ALLOCATABLE, DIMENSION(:) :: DFSQ
      CONTAINS
        SUBROUTINE ALLOC_DFACTSQ()
          USE DIM_MOD
          IF (ALLOCATED(DFSQ)) THEN
            DEALLOCATE(DFSQ)
          ENDIF
          ALLOCATE(DFSQ(0:2*NL_M-2))
        END SUBROUTINE ALLOC_DFACTSQ
      END MODULE DFACTSQ_MOD

C=======================================================================
      MODULE EIGEN_MOD
        IMPLICIT NONE
        INTEGER :: NE_EIG
        REAL :: E0_EIG
        REAL :: EFIN_EIG
        INTEGER :: I_VIB
        INTEGER :: I_MFP
      CONTAINS
        SUBROUTINE ALLOC_EIGEN()
          USE DIM_MOD
        END SUBROUTINE ALLOC_EIGEN
      END MODULE EIGEN_MOD

C=======================================================================
      MODULE FDIF_MOD
        IMPLICIT NONE
        REAL :: R1
        REAL :: R2
      CONTAINS
        SUBROUTINE ALLOC_FDIF()
          USE DIM_MOD
        END SUBROUTINE ALLOC_FDIF
      END MODULE FDIF_MOD

C=======================================================================
      MODULE FIXSCAN_MOD
        IMPLICIT NONE
        INTEGER :: N_FIXED
        INTEGER :: N_SCAN
        INTEGER :: IPH_1
        REAL :: FIX0
        REAL :: FIX1
        REAL :: SCAN0
        REAL :: SCAN1
      CONTAINS
        SUBROUTINE ALLOC_FIXSCAN()
          USE DIM_MOD
        END SUBROUTINE ALLOC_FIXSCAN
      END MODULE FIXSCAN_MOD

C=======================================================================
      MODULE FIXSCAN_A_MOD
        IMPLICIT NONE
        INTEGER :: N_FIXED_A
        INTEGER :: N_SCAN_A
        INTEGER :: IPH_1_A
        REAL :: FIX0_A
        REAL :: FIX1_A
        REAL :: SCAN0_A
        REAL :: SCAN1_A
      CONTAINS
        SUBROUTINE ALLOC_FIXSCAN_A()
          USE DIM_MOD
        END SUBROUTINE ALLOC_FIXSCAN_A
      END MODULE FIXSCAN_A_MOD

C=======================================================================
      MODULE LINLBD_MOD
        IMPLICIT NONE
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LBD
        INTEGER :: LBDMAX
        INTEGER, ALLOCATABLE, DIMENSION(:) :: NUMAX
      CONTAINS
        SUBROUTINE ALLOC_LINLBD()
          USE DIM_MOD
          IF (ALLOCATED(LBD)) THEN
            DEALLOCATE(LBD)
          ENDIF
          ALLOCATE(LBD(-N_MU_M:N_MU_M,0:N_NU_M))
          IF (ALLOCATED(NUMAX)) THEN
            DEALLOCATE(NUMAX)
          ENDIF
          ALLOCATE(NUMAX(NATP_M))
        END SUBROUTINE ALLOC_LINLBD
      END MODULE LINLBD_MOD

C=======================================================================
      MODULE MOYEN_A_MOD
        IMPLICIT NONE
        INTEGER :: IMOY_A
        INTEGER :: NDIR_A
        REAL :: ACCEPT_A
        INTEGER :: ICHKDIR_A
      CONTAINS
        SUBROUTINE ALLOC_MOYEN_A()
          USE DIM_MOD
        END SUBROUTINE ALLOC_MOYEN_A
      END MODULE MOYEN_A_MOD

C=======================================================================
      MODULE OUTFILES_MOD
        IMPLICIT NONE
        CHARACTER*24 :: OUTFILE1
        CHARACTER*24 :: OUTFILE2
        CHARACTER*24 :: OUTFILE3
        CHARACTER*24 :: OUTFILE4
      CONTAINS
        SUBROUTINE ALLOC_OUTFILES()
          USE DIM_MOD
        END SUBROUTINE ALLOC_OUTFILES
      END MODULE OUTFILES_MOD

C=======================================================================
      MODULE RA_MOD
        IMPLICIT NONE
        INTEGER :: I_NO
        INTEGER :: I_RA
        INTEGER, ALLOCATABLE, DIMENSION(:) :: N_RA
      CONTAINS
        SUBROUTINE ALLOC_RA()
          USE DIM_MOD
          IF (ALLOCATED(N_RA)) THEN
            DEALLOCATE(N_RA)
          ENDIF
          ALLOCATE(N_RA(NATP_M))
        END SUBROUTINE ALLOC_RA
      END MODULE RA_MOD

C=======================================================================
      MODULE SPECTRUM_MOD
        IMPLICIT NONE
        INTEGER, ALLOCATABLE, DIMENSION(:) :: I_SPECTRUM
      CONTAINS
        SUBROUTINE ALLOC_SPECTRUM()
          USE DIM_MOD
          IF (ALLOCATED(I_SPECTRUM)) THEN
            DEALLOCATE(I_SPECTRUM)
          ENDIF
          ALLOCATE(I_SPECTRUM(NE_M))
        END SUBROUTINE ALLOC_SPECTRUM
      END MODULE SPECTRUM_MOD

C=======================================================================
      MODULE DIRECT_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:,:) :: DIRANA
        REAL, ALLOCATABLE, DIMENSION(:,:) :: ANADIR
        REAL :: RTHEXT
        REAL :: RPHI
        REAL, ALLOCATABLE, DIMENSION(:) :: THETAR
        REAL, ALLOCATABLE, DIMENSION(:) :: PHIR
      CONTAINS
        SUBROUTINE ALLOC_DIRECT()
          USE DIM_MOD
          IF (ALLOCATED(DIRANA)) THEN
            DEALLOCATE(DIRANA)
          ENDIF
          ALLOCATE(DIRANA(3,49))
          IF (ALLOCATED(ANADIR)) THEN
            DEALLOCATE(ANADIR)
          ENDIF
          ALLOCATE(ANADIR(3,49))
          IF (ALLOCATED(THETAR)) THEN
            DEALLOCATE(THETAR)
          ENDIF
          ALLOCATE(THETAR(49))
          IF (ALLOCATED(PHIR)) THEN
            DEALLOCATE(PHIR)
          ENDIF
          ALLOCATE(PHIR(49))
        END SUBROUTINE ALLOC_DIRECT
      END MODULE DIRECT_MOD

C=======================================================================
      MODULE DIRECT_A_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:,:) :: DIRANA_A
        REAL, ALLOCATABLE, DIMENSION(:,:) :: ANADIR_A
        REAL :: RTHEXT_A
        REAL :: RPHI_A
        REAL, ALLOCATABLE, DIMENSION(:) :: THETAR_A
        REAL, ALLOCATABLE, DIMENSION(:) :: PHIR_A
      CONTAINS
        SUBROUTINE ALLOC_DIRECT_A()
          USE DIM_MOD
          IF (ALLOCATED(DIRANA_A)) THEN
            DEALLOCATE(DIRANA_A)
          ENDIF
          ALLOCATE(DIRANA_A(3,49))
          IF (ALLOCATED(ANADIR_A)) THEN
            DEALLOCATE(ANADIR_A)
          ENDIF
          ALLOCATE(ANADIR_A(3,49))
          IF (ALLOCATED(THETAR_A)) THEN
            DEALLOCATE(THETAR_A)
          ENDIF
          ALLOCATE(THETAR_A(49))
          IF (ALLOCATED(PHIR_A)) THEN
            DEALLOCATE(PHIR_A)
          ENDIF
          ALLOCATE(PHIR_A(49))
        END SUBROUTINE ALLOC_DIRECT_A
      END MODULE DIRECT_A_MOD

C=======================================================================
      MODULE PATH_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:) :: NPATH
        REAL, ALLOCATABLE, DIMENSION(:) :: NPATH2
        INTEGER :: NTHOF
        REAL, ALLOCATABLE, DIMENSION(:) :: NPMA
        REAL, ALLOCATABLE, DIMENSION(:) :: NPMI
      CONTAINS
        SUBROUTINE ALLOC_PATH()
          USE DIM_MOD
          IF (ALLOCATED(NPATH)) THEN
            DEALLOCATE(NPATH)
          ENDIF
          ALLOCATE(NPATH(0:NDIF_M))
          IF (ALLOCATED(NPATH2)) THEN
            DEALLOCATE(NPATH2)
          ENDIF
          ALLOCATE(NPATH2(0:NDIF_M))
          IF (ALLOCATED(NPMA)) THEN
            DEALLOCATE(NPMA)
          ENDIF
          ALLOCATE(NPMA(0:NDIF_M))
          IF (ALLOCATED(NPMI)) THEN
            DEALLOCATE(NPMI)
          ENDIF
          ALLOCATE(NPMI(0:NDIF_M))
        END SUBROUTINE ALLOC_PATH
      END MODULE PATH_MOD

C=======================================================================
      MODULE ROT_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: RLM01
      CONTAINS
        SUBROUTINE ALLOC_ROT()
          USE DIM_MOD
          IF (ALLOCATED(RLM01)) THEN
            DEALLOCATE(RLM01)
          ENDIF
          ALLOCATE(RLM01(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1))
        END SUBROUTINE ALLOC_ROT
      END MODULE ROT_MOD

C=======================================================================
      MODULE TESTPA_MOD
        IMPLICIT NONE
        INTEGER, ALLOCATABLE, DIMENSION(:) :: IT
        INTEGER, ALLOCATABLE, DIMENSION(:) :: IN
        INTEGER :: IJ
      CONTAINS
        SUBROUTINE ALLOC_TESTPA()
          USE DIM_MOD
          IF (ALLOCATED(IT)) THEN
            DEALLOCATE(IT)
          ENDIF
          ALLOCATE(IT(0:NDIF_M))
          IF (ALLOCATED(IN)) THEN
            DEALLOCATE(IN)
          ENDIF
          ALLOCATE(IN(0:NDIF_M))
        END SUBROUTINE ALLOC_TESTPA
      END MODULE TESTPA_MOD

C=======================================================================
      MODULE TESTPB_MOD
        IMPLICIT NONE
        COMPLEX :: RHO01
        REAL :: TH01
        REAL :: PHI01
      CONTAINS
        SUBROUTINE ALLOC_TESTPB()
          USE DIM_MOD
        END SUBROUTINE ALLOC_TESTPB
      END MODULE TESTPB_MOD

C=======================================================================
      MODULE TLDW_MOD
        IMPLICIT NONE
        COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: TLT
      CONTAINS
        SUBROUTINE ALLOC_TLDW()
          USE DIM_MOD
          IF (ALLOCATED(TLT)) THEN
            DEALLOCATE(TLT)
          ENDIF
          ALLOCATE(TLT(0:NT_M,4,NATM,NE_M))
        END SUBROUTINE ALLOC_TLDW
      END MODULE TLDW_MOD

C=======================================================================
      MODULE VARIA_MOD
        IMPLICIT NONE
        COMPLEX, ALLOCATABLE, DIMENSION(:) :: CEX
        COMPLEX, ALLOCATABLE, DIMENSION(:) :: CEXDW
        REAL, ALLOCATABLE, DIMENSION(:) :: DW
      CONTAINS
        SUBROUTINE ALLOC_VARIA()
          USE DIM_MOD
          IF (ALLOCATED(CEX)) THEN
            DEALLOCATE(CEX)
          ENDIF
          ALLOCATE(CEX(0:NDIF_M))
          IF (ALLOCATED(CEXDW)) THEN
            DEALLOCATE(CEXDW)
          ENDIF
          ALLOCATE(CEXDW(0:NDIF_M))
          IF (ALLOCATED(DW)) THEN
            DEALLOCATE(DW)
          ENDIF
          ALLOCATE(DW(0:NDIF_M))
        END SUBROUTINE ALLOC_VARIA
      END MODULE VARIA_MOD

C=======================================================================
      MODULE LBD_MOD
        IMPLICIT NONE
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LBDM
      CONTAINS
        SUBROUTINE ALLOC_LBD()
          USE DIM_MOD
          INTEGER I, J
          IF (ALLOCATED(LBDM)) THEN
            DEALLOCATE(LBDM)
          ENDIF
          ALLOCATE(LBDM(2,NDIF_M))
          DO I=1,NDIF_M
            DO J=1,2
              LBDM(J,I) = 0
            ENDDO
          ENDDO
        END SUBROUTINE ALLOC_LBD
      END MODULE LBD_MOD

C=======================================================================
      MODULE SCATMAT_MOD
        IMPLICIT NONE
        COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: F21
      CONTAINS
        SUBROUTINE ALLOC_SCATMAT()
          USE DIM_MOD
          IF (ALLOCATED(F21)) THEN
            DEALLOCATE(F21)
          ENDIF
          ALLOCATE(F21(NSPIN2_M,NLAMBDA_M,NLAMBDA_M,NDIF_M))
        END SUBROUTINE ALLOC_SCATMAT
      END MODULE SCATMAT_MOD

C=======================================================================
      MODULE EXTREM_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:) :: FMIN
        REAL, ALLOCATABLE, DIMENSION(:) :: FMAX
        INTEGER :: IREF
      CONTAINS
        SUBROUTINE ALLOC_EXTREM()
          USE DIM_MOD
          IF (ALLOCATED(FMIN)) THEN
            DEALLOCATE(FMIN)
          ENDIF
          ALLOCATE(FMIN(0:NDIF_M))
          IF (ALLOCATED(FMAX)) THEN
            DEALLOCATE(FMAX)
          ENDIF
          ALLOCATE(FMAX(0:NDIF_M))
        END SUBROUTINE ALLOC_EXTREM
      END MODULE EXTREM_MOD

C=======================================================================
      MODULE PRINTP_MOD
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:) :: FMN
        REAL, ALLOCATABLE, DIMENSION(:) :: PATH
        REAL, ALLOCATABLE, DIMENSION(:) :: DMN
        INTEGER, ALLOCATABLE, DIMENSION(:) :: JON
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: JPON
      CONTAINS
        SUBROUTINE ALLOC_PRINTP()
          USE DIM_MOD
          IF (ALLOCATED(FMN)) THEN
            DEALLOCATE(FMN)
          ENDIF
          ALLOCATE(FMN(NPATH_M))
          IF (ALLOCATED(PATH)) THEN
            DEALLOCATE(PATH)
          ENDIF
          ALLOCATE(PATH(NPATH_M))
          IF (ALLOCATED(DMN)) THEN
            DEALLOCATE(DMN)
          ENDIF
          ALLOCATE(DMN(NPATH_M))
          IF (ALLOCATED(JON)) THEN
            DEALLOCATE(JON)
          ENDIF
          ALLOCATE(JON(NPATH_M))
          IF (ALLOCATED(JPON)) THEN
            DEALLOCATE(JPON)
          ENDIF
          ALLOCATE(JPON(NPATH_M,NDIF_M))
        END SUBROUTINE ALLOC_PRINTP
      END MODULE PRINTP_MOD
C=======================================================================
      MODULE ROT_CUB_MOD
        IMPLICIT NONE
        REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: R_PIS2
      CONTAINS
        SUBROUTINE ALLOC_ROT_CUB()
          USE DIM_MOD
          IF (ALLOCATED(R_PIS2)) THEN
            DEALLOCATE(R_PIS2)
          ENDIF
          ALLOCATE(R_PIS2(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1))
        END SUBROUTINE ALLOC_ROT_CUB
      END MODULE ROT_CUB_MOD
C=======================================================================
      MODULE SYMMOP_MOD
        IMPLICIT NONE
        REAL*4, ALLOCATABLE, DIMENSION(:) :: ZL
        COMPLEX, ALLOCATABLE, DIMENSION(:) :: ZM1
        COMPLEX, ALLOCATABLE, DIMENSION(:) :: ZM2
        INTEGER, ALLOCATABLE, DIMENSION(:) :: IZ
      CONTAINS
        SUBROUTINE ALLOC_SYMMOP()
          USE DIM_MOD
          IF (ALLOCATED(ZL)) THEN
            DEALLOCATE(ZL)
          ENDIF
          ALLOCATE(ZL(64))
          IF (ALLOCATED(ZM1)) THEN
            DEALLOCATE(ZM1)
          ENDIF
          ALLOCATE(ZM1(64))
          IF (ALLOCATED(ZM2)) THEN
            DEALLOCATE(ZM2)
          ENDIF
          ALLOCATE(ZM2(64))
          IF (ALLOCATED(IZ)) THEN
            DEALLOCATE(IZ)
          ENDIF
          ALLOCATE(IZ(64))
        END SUBROUTINE ALLOC_SYMMOP
      END MODULE SYMMOP_MOD
C=======================================================================
      MODULE TAU_PROT_MOD
        IMPLICIT NONE
        REAL*4, ALLOCATABLE, DIMENSION(:) :: Z_L_P
        COMPLEX, ALLOCATABLE, DIMENSION(:) :: Z_M_P
        INTEGER, ALLOCATABLE, DIMENSION(:) :: I_Z_P
        INTEGER, ALLOCATABLE, DIMENSION(:) :: ISTEP_L
        INTEGER, ALLOCATABLE, DIMENSION(:) :: ISTEP_M
        INTEGER, ALLOCATABLE, DIMENSION(:) :: I_LM
        INTEGER, ALLOCATABLE, DIMENSION(:) :: I_REL_MP
      CONTAINS
        SUBROUTINE ALLOC_TAU_PROT()
          USE DIM_MOD
          IF (ALLOCATED(Z_L_P)) THEN
            DEALLOCATE(Z_L_P)
          ENDIF
          ALLOCATE(Z_L_P(NATP_M))
          IF (ALLOCATED(Z_M_P)) THEN
            DEALLOCATE(Z_M_P)
          ENDIF
          ALLOCATE(Z_M_P(NATP_M))
          IF (ALLOCATED(I_Z_P)) THEN
            DEALLOCATE(I_Z_P)
          ENDIF
          ALLOCATE(I_Z_P(NATP_M))
          IF (ALLOCATED(ISTEP_L)) THEN
            DEALLOCATE(ISTEP_L)
          ENDIF
          ALLOCATE(ISTEP_L(NATP_M))
          IF (ALLOCATED(ISTEP_M)) THEN
            DEALLOCATE(ISTEP_M)
          ENDIF
          ALLOCATE(ISTEP_M(NATP_M))
          IF (ALLOCATED(I_LM)) THEN
            DEALLOCATE(I_LM)
          ENDIF
          ALLOCATE(I_LM(NATP_M))
          IF (ALLOCATED(I_REL_MP)) THEN
            DEALLOCATE(I_REL_MP)
          ENDIF
          ALLOCATE(I_REL_MP(NATP_M))
        END SUBROUTINE ALLOC_TAU_PROT
      END MODULE TAU_PROT_MOD
C=======================================================================
      MODULE TAUSYMMOP_MOD
        IMPLICIT NONE
        REAL*4, ALLOCATABLE, DIMENSION(:,:) :: Z_L
        COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: Z_M1
        COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: Z_M2
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: I_Z
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISYM
        INTEGER, ALLOCATABLE, DIMENSION(:) :: NSYM_G
        INTEGER :: NSIZE_GR
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: S_INV
      CONTAINS
        SUBROUTINE ALLOC_TAUSYMMOP()
          USE DIM_MOD
          IF (ALLOCATED(Z_L)) THEN
            DEALLOCATE(Z_L)
          ENDIF
          ALLOCATE(Z_L(NAT_EQ_M,NATP_M))
          IF (ALLOCATED(Z_M1)) THEN
            DEALLOCATE(Z_M1)
          ENDIF
          ALLOCATE(Z_M1(NAT_EQ_M,NATP_M))
          IF (ALLOCATED(Z_M2)) THEN
            DEALLOCATE(Z_M2)
          ENDIF
          ALLOCATE(Z_M2(NAT_EQ_M,NATP_M))
          IF (ALLOCATED(I_Z)) THEN
            DEALLOCATE(I_Z)
          ENDIF
          ALLOCATE(I_Z(NAT_EQ_M,NATP_M))
          IF (ALLOCATED(ISYM)) THEN
            DEALLOCATE(ISYM)
          ENDIF
          ALLOCATE(ISYM(NAT_EQ_M,NATP_M))
          IF (ALLOCATED(NSYM_G)) THEN
            DEALLOCATE(NSYM_G)
          ENDIF
          ALLOCATE(NSYM_G(48))
          IF (ALLOCATED(S_INV)) THEN
            DEALLOCATE(S_INV)
          ENDIF
          ALLOCATE(S_INV(64,3,3))
        END SUBROUTINE ALLOC_TAUSYMMOP
      END MODULE TAUSYMMOP_MOD
C=======================================================================
      MODULE RENORM_MOD
        IMPLICIT NONE
        INTEGER :: I_REN, N_REN
        REAL    :: REN_R, REN_I
      CONTAINS
        SUBROUTINE ALLOC_RENORM()
          USE DIM_MOD
        END SUBROUTINE ALLOC_RENORM
      END MODULE RENORM_MOD
C=======================================================================
      MODULE C_RENORM_MOD
        IMPLICIT NONE
        COMPLEX, ALLOCATABLE, DIMENSION(:) :: C_REN
      CONTAINS
        SUBROUTINE ALLOC_C_RENORM()
          USE DIM_MOD
          IF (ALLOCATED(C_REN)) THEN
            DEALLOCATE(C_REN)
          ENDIF
          ALLOCATE(C_REN(0:NDIF_M))
        END SUBROUTINE ALLOC_C_RENORM
      END MODULE C_RENORM_MOD
