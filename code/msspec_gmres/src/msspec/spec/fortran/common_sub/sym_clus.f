C
C=======================================================================
C
      SUBROUTINE SYM_CLUS(CUNIT,COORD,CHEM_OLD,NZAT,NTYP,IPHA)
C
C  This subroutine reorganizes the cluster into equivalence classes
C    taking advantage of its symmetry properties. The corresponding
C    symmetry operations are stored for further use.
C
C                           Written by Mihai Gavaza      : 15 Feb 2000
C                           Modified by Didier Sebilleau : 05 Nov 2003
C
C      INCLUDE 'spec.inc'
      USE DIM_MOD
C
      USE ATOMS_MOD, NZ_AT => NZAT
      USE CLUSLIM_MOD, Z_PLAN => VAL
      USE COOR_MOD
      USE OUTFILES_MOD
      USE OUTUNITS_MOD
      USE ROT_CUB_MOD
      USE SYMMOP_MOD
      USE TAU_PROT_MOD
      USE TAUSYMMOP_MOD
      USE TESTS_MOD
C
C
C
      COMPLEX ONEC,IC,JC1,JC2,JC4,JC5
      COMPLEX ONEC_,IC_,JC1_,JC2_,JC4_,JC5_
C
      PARAMETER(ONEC=(1.,0.),IC=(0.,1.),JC1=(0.5,0.866025))
      PARAMETER(ONEC_=(-1.,0.),IC_=(0.,-1.),JC1_=(-0.5,-0.866025))
      PARAMETER(JC2=JC1*JC1,JC4=JC2*JC2,JC5=JC4*JC1)
      PARAMETER(JC2_=ONEC_*JC2,JC4_=ONEC_*JC4,JC5_=ONEC_*JC5)
C
      REAL*4 COORD(3,NATCLU_M),COORD1(3,NATCLU_M),SYM_AT1(3,NATCLU_M)
      REAL*4 DNBR(NATCLU_M,NATCLU_M),DIST0(NATCLU_M),DIST1(NATCLU_M)
      REAL*4 S_M(3,3,64),S_MUL(3,3)
C
      INTEGER NSYM_M(64),CONT_G(48,32),SIZE_G(32)
      INTEGER INV_P(NATP_M,64),INV_PT(NATP_M,64),NAT_SYM(4)
      INTEGER I_SET(NATP_M),GR(NATP_M)
      INTEGER NEQAT(NATCLU_M),NBRZ(NATCLU_M,NATCLU_M),NZAT(NATCLU_M)
      INTEGER IUSED(NATCLU_M),NTYP(NATCLU_M),NGAIN_B(NATP_M)
      INTEGER NQAT(NATCLU_M),I_33(32),IS_WRONG(64),OLD_G(32),NEW_G(32)
      INTEGER N_NEW_OLD(NATCLU_M)
C
      CHARACTER*2 CHEM_OLD(NATCLU_M)
      CHARACTER*3 NAME_G(0:32)
      CHARACTER*4 ANG_ROT(4)
      CHARACTER*5 NAME_S(64)
      CHARACTER*8 AT_ADD(2)
C
C
      LOGICAL EQUIV,MATCH
C
C  Matrices of the 64 symmetry operations in the form S_M(I,J)
C     I is chosen as the line index and J as the column one
C
C  Matrices for the rotations of the group Oh
C
C     The greek indices for the cubic rotations have been written as :
C
C                     alpha ----> l
C                     beta  ----> m
C                     gamma ----> n
C                     delta ----> o
C
      DATA ((S_M(I,J,1), J=1,3), I=1,3)  /1.,0.,0.,0.,1.,0.,0.,0.,1./  
     &        ! E
      DATA ((S_M(I,J,2), J=1,3), I=1,3)  /1.,0.,0.,0.,-1.,0.,0.,0.,-1./
     &        ! C2X
      DATA ((S_M(I,J,3), J=1,3), I=1,3)  /-1.,0.,0.,0.,1.,0.,0.,0.,-1./
     &        ! C2Y
      DATA ((S_M(I,J,4), J=1,3), I=1,3)  /-1.,0.,0.,0.,-1.,0.,0.,0.,1./
     &        ! C2Z
      DATA ((S_M(I,J,5), J=1,3), I=1,3)  /0.,0.,1.,1.,0.,0.,0.,1.,0./  
     &        ! C3l
      DATA ((S_M(I,J,6), J=1,3), I=1,3)  /0.,-1.,0.,0.,0.,1.,-1.,0.,0./
     &        ! C3m
      DATA ((S_M(I,J,7), J=1,3), I=1,3)  /0.,0.,-1.,1.,0.,0.,0.,-1.,0./
     &        ! C3n
      DATA ((S_M(I,J,8), J=1,3), I=1,3)  /0.,-1.,0.,0.,0.,-1.,1.,0.,0./
     &        ! C3o
      DATA ((S_M(I,J,9), J=1,3), I=1,3)  /0.,1.,0.,0.,0.,1.,1.,0.,0./  
     &        ! C3l2
      DATA ((S_M(I,J,10), J=1,3), I=1,3) /0.,0.,-1.,-1.,0.,0.,0.,1.,0./
     &        ! C3m2
      DATA ((S_M(I,J,11), J=1,3), I=1,3) /0.,1.,0.,0.,0.,-1.,-1.,0.,0./
     &        ! C3n2
      DATA ((S_M(I,J,12), J=1,3), I=1,3) /0.,0.,1.,-1.,0.,0.,0.,-1.,0./
     &        ! C3o2
      DATA ((S_M(I,J,13), J=1,3), I=1,3) /1.,0.,0.,0.,0.,-1.,0.,1.,0./ 
     &        ! C4X
      DATA ((S_M(I,J,14), J=1,3), I=1,3) /0.,0.,1.,0.,1.,0.,-1.,0.,0./ 
     &        ! C4Y
      DATA ((S_M(I,J,15), J=1,3), I=1,3) /0.,-1.,0.,1.,0.,0.,0.,0.,1./ 
     &        ! C4Z
      DATA ((S_M(I,J,16), J=1,3), I=1,3) /1.,0.,0.,0.,0.,1.,0.,-1.,0./ 
     &        ! C4X3
      DATA ((S_M(I,J,17), J=1,3), I=1,3) /0.,0.,-1.,0.,1.,0.,1.,0.,0./ 
     &        ! C4Y3
      DATA ((S_M(I,J,18), J=1,3), I=1,3) /0.,1.,0.,-1.,0.,0.,0.,0.,1./ 
     &        ! C4Z3
      DATA ((S_M(I,J,19), J=1,3), I=1,3) /0.,1.,0.,1.,0.,0.,0.,0.,-1./ 
     &        ! C2a
      DATA ((S_M(I,J,20), J=1,3), I=1,3) /0.,-1.,0.,-1.,0.,0.,0.,0.,-1.
     &/       ! C2b
      DATA ((S_M(I,J,21), J=1,3), I=1,3) /0.,0.,1.,0.,-1.,0.,1.,0.,0./ 
     &        ! C2c
      DATA ((S_M(I,J,22), J=1,3), I=1,3) /0.,0.,-1.,0.,-1.,0.,-1.,0.,0.
     &/       ! C2d
      DATA ((S_M(I,J,23), J=1,3), I=1,3) /-1.,0.,0.,0.,0.,-1.,0.,-1.,0.
     &/       ! C2e
      DATA ((S_M(I,J,24), J=1,3), I=1,3) /-1.,0.,0.,0.,0.,1.,0.,1.,0./ 
     &        ! C2f
C
C  Matrices for the rotations of the group D6h
C
      DATA ((S_M(I,J,25), J=1,3), I=1,3) /-0.5,-0.866025,0.,0.866025,-
     &0.5,0.,0.,0.,1./   ! C3Z
      DATA ((S_M(I,J,26), J=1,3), I=1,3) /-0.5,0.866025,0.,-0.866025,-
     &0.5,0.,0.,0.,1./   ! C3Z2
      DATA ((S_M(I,J,27), J=1,3), I=1,3) /0.5,-0.866025,0.,0.866025,0.
     &5,0.,0.,0.,1./     ! C6Z
      DATA ((S_M(I,J,28), J=1,3), I=1,3) /0.5,0.866025,0.,-0.866025,0.
     &5,0.,0.,0.,1./     ! C6Z5
      DATA ((S_M(I,J,29), J=1,3), I=1,3) /-0.5,-0.866025,0.,-0.866025,
     &0.5,0.,0.,0.,-1./  ! C2A
      DATA ((S_M(I,J,30), J=1,3), I=1,3) /-0.5,0.866025,0.,0.866025,0.
     &5,0.,0.,0.,-1./    ! C2B
      DATA ((S_M(I,J,31), J=1,3), I=1,3) /0.5,-0.866025,0.,-0.866025,-
     &0.5,0.,0.,0.,-1./  ! C2C
      DATA ((S_M(I,J,32), J=1,3), I=1,3) /0.5,0.866025,0.,0.866025,-0.
     &5,0.,0.,0.,-1./    ! C2D
C
C  Matrices for the roto-inversions of the group Oh
C
      DATA ((S_M(I,J,33), J=1,3), I=1,3) /-1.,0.,0.,0.,-1.,0.,0.,0.,-1.
     &/       ! I
      DATA ((S_M(I,J,34), J=1,3), I=1,3)  /-1.,0.,0.,0.,1.,0.,0.,0.,1./
     &         ! IC2X
      DATA ((S_M(I,J,35), J=1,3), I=1,3)  /1.,0.,0.,0.,-1.,0.,0.,0.,1./
     &         ! IC2Y
      DATA ((S_M(I,J,36), J=1,3), I=1,3)  /1.,0.,0.,0.,1.,0.,0.,0.,-1./
     &         ! IC2Z
      DATA ((S_M(I,J,37), J=1,3), I=1,3)  /0.,0.,-1.,-1.,0.,0.,0.,-1.,
     &0./       ! IC3l
      DATA ((S_M(I,J,38), J=1,3), I=1,3)/0.,1.,0.,0.,0.,-1.,1.,0.,0./  
     &       ! IC3m
      DATA ((S_M(I,J,39), J=1,3), I=1,3)  /0.,0.,1.,-1.,0.,0.,0.,1.,0./
     &         ! IC3n
      DATA ((S_M(I,J,40), J=1,3), I=1,3) /0.,1.,0.,0.,0.,1.,-1.,0.,0./ 
     &        ! IC3o
      DATA ((S_M(I,J,41), J=1,3), I=1,3) /0.,-1.,0.,0.,0.,-1.,-1.,0.,0.
     &/       ! IC3l2
      DATA ((S_M(I,J,42), J=1,3), I=1,3) /0.,0.,1.,1.,0.,0.,0.,-1.,0./ 
     &        ! IC3m2
      DATA ((S_M(I,J,43), J=1,3), I=1,3) /0.,-1.,0.,0.,0.,1.,1.,0.,0./ 
     &        ! IC3n2
      DATA ((S_M(I,J,44), J=1,3), I=1,3) /0.,0.,-1.,1.,0.,0.,0.,1.,0./ 
     &        ! IC3o2
      DATA ((S_M(I,J,45), J=1,3), I=1,3) /-1.,0.,0.,0.,0.,1.,0.,-1.,0./
     &        ! IC4X
      DATA ((S_M(I,J,46), J=1,3), I=1,3) /0.,0.,-1.,0.,-1.,0.,1.,0.,0./
     &        ! IC4Y
      DATA ((S_M(I,J,47), J=1,3), I=1,3) /0.,1.,0.,-1.,0.,0.,0.,0.,-1./
     &        ! IC4Z
      DATA ((S_M(I,J,48), J=1,3), I=1,3) /-1.,0.,0.,0.,0.,-1.,0.,1.,0./
     &        ! IC4X3
      DATA ((S_M(I,J,49), J=1,3), I=1,3) /0.,0.,1.,0.,-1.,0.,-1.,0.,0./
     &        ! IC4Y3
      DATA ((S_M(I,J,50), J=1,3), I=1,3) /0.,-1.,0.,1.,0.,0.,0.,0.,-1./
     &        ! IC4Z3
      DATA ((S_M(I,J,51), J=1,3), I=1,3) /0.,-1.,0.,-1.,0.,0.,0.,0.,1./
     &        ! IC2a
      DATA ((S_M(I,J,52), J=1,3), I=1,3) /0.,1.,0.,1.,0.,0.,0.,0.,1./  
     &        ! IC2b
      DATA ((S_M(I,J,53), J=1,3), I=1,3) /0.,0.,-1.,0.,1.,0.,-1.,0.,0./
     &        ! IC2c
      DATA ((S_M(I,J,54), J=1,3), I=1,3) /0.,0.,1.,0.,1.,0.,1.,0.,0./  
     &        ! IC2d
      DATA ((S_M(I,J,55), J=1,3), I=1,3) /1.,0.,0.,0.,0.,1.,0.,1.,0./  
     &        ! IC2e
      DATA ((S_M(I,J,56), J=1,3), I=1,3) /1.,0.,0.,0.,0.,-1.,0.,-1.,0./
     &        ! IC2f
C
C  Matrices for the roto-inversions of the group D6h
C
      DATA ((S_M(I,J,57), J=1,3), I=1,3) /0.5,0.866025,0.,-0.866025,0.
     &5,0.,0.,0.,-1./    ! IC3Z
      DATA ((S_M(I,J,58), J=1,3), I=1,3) /0.5,-0.866025,0.,0.866025,0.
     &5,0.,0.,0.,-1./    ! IC3Z2
      DATA ((S_M(I,J,59), J=1,3), I=1,3) /-0.5,0.866025,0.,-0.866025,-
     &0.5,0.,0.,0.,-1./  ! IC6Z
      DATA ((S_M(I,J,60), J=1,3), I=1,3) /-0.5,-0.866025,0.,0.866025,-
     &0.5,0.,0.,0.,-1./  ! IC6Z5
      DATA ((S_M(I,J,61), J=1,3), I=1,3) /0.5,0.866025,0.,0.866025,-0.
     &5,0.,0.,0.,1./     ! IC2A
      DATA ((S_M(I,J,62), J=1,3), I=1,3)/0.5,-0.866025,0.,-0.866025,-0.
     &5,0.,0.,0.,1./   ! IC2B
      DATA ((S_M(I,J,63), J=1,3), I=1,3) /-0.5,0.866025,0.,0.866025,0.
     &5,0.,0.,0.,1./     ! IC2C
      DATA ((S_M(I,J,64), J=1,3), I=1,3) /-0.5,-0.866025,0.,-0.866025,
     &0.5,0.,0.,0.,1./   ! IC2D
C
C  For a given symmetry operation S, IZ can have three values :
C
C         IZ =  1  ---->  delta_{L, L'} in S_{L L'}
C         IZ =  0  ---->  delta_{l, l'} only in S_{L L'}
C         IZ = -1  ---->  delta_{L,-L'} in S_{L L'}
C
      IZ = (/1,-1,-1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,-1,-1,0,0,0,0,1,1,
     &1,1,-1,-1,-1,-1,1,-1,-1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,-1,-1,0,0,
     &0,0,1,1,1,1,-1,-1,-1,-1/)
C
C  For a given symmetry operation S, S_{L L'} is proportional to
C                      ZL**L and ZM1**M (and ZM2**M' if IZ=0)
C
      ZL = (/1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
     &-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,-1.,-1.,-1.,-1.,-1.,1.,1.,-1.,-1.
     &,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,1.,1.,-1.,-
     &1.,-1.,-1.,-1.,-1.,-1.,-1.,1.,1.,1.,1./)
C
      ZM1 = (/ONEC,ONEC,ONEC_,ONEC_,ONEC,IC_,ONEC_,IC,IC_,ONEC_,IC,
     &ONEC,IC,ONEC,IC_,IC_,ONEC_,IC,IC_,IC,ONEC,ONEC_,IC,IC_,JC4,JC2,
     &JC5,JC1,JC5_,JC1_,JC4_,JC2_,ONEC,ONEC,ONEC_,ONEC_,ONEC,IC_,ONEC_,
     &IC,IC_,ONEC_,IC,ONEC,IC,ONEC,IC_,IC_,ONEC_,IC,IC_,IC,ONEC,ONEC_,
     &IC,IC_,JC4,JC2,JC5,JC1,JC5_,JC1_,JC4_,JC2_/)
      ZM2 = (/ONEC,ONEC,ONEC,ONEC,IC_,ONEC,IC,ONEC_,ONEC_,IC_,ONEC,
     &IC,IC_,ONEC,ONEC,IC,ONEC_,ONEC,ONEC,ONEC,ONEC_,ONEC,IC,IC_,ONEC,
     &ONEC,ONEC,ONEC,ONEC,ONEC,ONEC,ONEC,ONEC,ONEC,ONEC,ONEC,IC_,ONEC,
     &IC,ONEC_,ONEC_,IC_,ONEC,IC,IC_,ONEC,ONEC,ONEC,ONEC,ONEC,ONEC,
     &ONEC,ONEC_,ONEC,IC,IC_,ONEC,ONEC,IC,ONEC_,ONEC,ONEC,ONEC,ONEC/)
C
C  Name of the crystallographic point-groups
C
      DATA NAME_G /'---',
     &             ' C1',' Ci',' C2','C1h','C2h',' D2','C2v','D2h',
     &             ' C4',' S4','C4h',' D4','C4v','D2d','D4h',' C3',
     &             'C3i',' D3','C3v','D3d',' C6','C3h','C6h',' D6',
     &             'C6v','D3h','D6h','  T',' Th',' Td','  O',' Oh'/
C
C  Content of the crystallographic point-groups : CONT_G(JROT,JGROUP)
C
C     In some cases, two contents are given. They correspond to a rotation of
C         the x and y axes about 0z by pi/4 in the cube (D2,C2v,D2h,D2d)
C         or pi/6 in the hexagon (D3,C3v,D3d,D3h). In this case, x and y
C          are respectively transformed into a and b in the first case,
C           and D and A in the second. The cube is invariant by any pi/2
C             rotation and the hexagon by any pi/3 rotation about 0z.
C
C
      DATA (CONT_G(I,1), I=1,48) /1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,     
     &    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0/            ! C1
      DATA (CONT_G(I,2), I=1,48) /1,33,0,0,0,0,0,0,0,0,0,0,0,0,0,0,    
     &    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0/            ! Ci
      DATA (CONT_G(I,3), I=1,48) /1,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/     
     &       ! C2
      DATA (CONT_G(I,4), I=1,48) /1,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/     
     &       ! C1h
      DATA (CONT_G(I,5), I=1,48) /1,4,33,36,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   
     &         ! C2h
      DATA (CONT_G(I,6), I=1,48) /1,2,3,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,4,19,20/   
     &       ! D2
      DATA (CONT_G(I,7), I=1,48) /1,4,34,35,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,4,51,52/ 
     &         ! C2v
      DATA (CONT_G(I,8), I=1,48) /1,2,3,4,33,34,35,36,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,4,19,20,33,36,
     &51,52/      ! D2h
      DATA (CONT_G(I,9), I=1,48) /1,4,15,18,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   
     &         ! C4
      DATA (CONT_G(I,10), I=1,48)/1,4,47,50,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   
     &         ! S4
      DATA (CONT_G(I,11), I=1,48)/1,4,15,18,33,36,47,50,0,0,0,0,0,0,
     *                     0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &                     0/            ! C4h
      DATA (CONT_G(I,12), I=1,48)/1,2,3,4,15,18,19,20,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/ 
     &           ! D4
      DATA (CONT_G(I,13), I=1,48)/1,4,15,18,34,35,51,52,0,0,0,0,0,0,
     *                     0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &                     0/            ! C4v
      DATA (CONT_G(I,14), I=1,48)/1,2,3,4,47,50,51,52,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,4,19,20,34,35,
     &47,50/      ! D2d
      DATA (CONT_G(I,15), I=1,48)/1,2,3,4,15,18,19,20,33,34,35,36,47,
     *                     50,51,52,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0/            ! D4h
      DATA (CONT_G(I,16), I=1,48)/1,25,26,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   
     &         ! C3
      DATA (CONT_G(I,17), I=1,48)/1,25,26,33,57,58,0,0,0,0,0,0,0,0,
     *                     0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &                     0/            ! C3i
      DATA (CONT_G(I,18), I=1,48)/1,3,25,26,31,32,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,25,26,29,
     &30/        ! D3
      DATA (CONT_G(I,19), I=1,48)/1,25,26,34,61,62,0,0,0,0,0,0,0,0,
     *                     0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,25,26,
     &                     35,63,64/       ! C3v
      DATA (CONT_G(I,20), I=1,48)/1,3,25,26,31,32,33,35,57,58,63,64,
     *                     0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,25,26,29,30,
     &                     33,34,57,58,61,62/  ! D3d
      DATA (CONT_G(I,21), I=1,48)/1,4,25,26,27,28,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/ 
     &           ! C6
      DATA (CONT_G(I,22), I=1,48)/1,25,26,36,59,60,0,0,0,0,0,0,0,0,
     *                     0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &                     0/            ! C3h
      DATA (CONT_G(I,23), I=1,48)/1,4,25,26,27,28,33,36,57,58,59,60,
     *                     0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0/            ! C6h
      DATA (CONT_G(I,24), I=1,48)/1,2,3,4,25,26,27,28,29,30,31,32,0,
     *                     0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0/            ! D6
      DATA (CONT_G(I,25), I=1,48)/1,4,25,26,27,28,34,35,61,62,63,64,
     *                     0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0/            ! C6v
      DATA (CONT_G(I,26), I=1,48)/1,3,25,26,31,32,34,36,59,60,61,62,
     *                     0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,25,26,29,30,
     &                     35,36,59,60,63,64/  ! D3h
      DATA (CONT_G(I,27), I=1,48)/1,2,3,4,25,26,27,28,29,30,31,32,33,
     *                     34,35,36,57,58,59,60,
     &                     61,62,63,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0,0/            ! D6h
      DATA (CONT_G(I,28), I=1,48)/1,2,3,4,5,6,7,8,9,10,11,12,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   
     &         ! T
      DATA (CONT_G(I,29), I=1,48)/1,2,3,4,5,6,7,8,9,10,11,12,33,34,
     *                     35,36,37,38,39,40,41,
     &                     42,43,44,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0/            ! Th
      DATA (CONT_G(I,30), I=1,48)/1,2,3,4,5,6,7,8,9,10,11,12,45,46,
     *                     47,48,49,50,51,52,53,
     &                     54,55,56,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0/            ! Td
      DATA (CONT_G(I,31), I=1,48)/1,2,3,4,5,6,7,8,9,10,11,12,13,14,
     *                     15,16,17,18,19,20,21,
     &                     22,23,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &                     0,0,0,0,0,0,0/            ! O
      DATA (CONT_G(I,32), I=1,48)/1,2,3,4,5,6,7,8,9,10,11,12,13,14,
     *                15,16,17,18,19,20,21,22,23,24,33,
     &                34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,
     &                50,51,52,53,54,55,56/ ! Oh
C
C  Size of the point-groups
C
      DATA SIZE_G /1,2,2,2,4,4,4,8,4,4,8,8,8,8,16,3,6,6,6,12,6,6,12,12,
     &12,12,24,12,24,24,24,48/
C
C  Groups containing the space inversion I
C
      DATA I_33   /0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,0,1,0,0,1,0,0,1,0,0,0,
     &1,0,1,0,0,1/
C
C  Name of the symmetry operations
C
C  Note : the three-fold rotation axes alpha, beta, gamma and delta
C         of the cubic groups have been coded respectively
C         l, m, n and o
C
      DATA NAME_S /'    E','  C2X','  C2Y','  C2Z','  C3l','  C3m',
     &             '  C3n','  C3o',' C3l2',' C3m2',' C3n2',' C3o2',
     &             '  C4X','  C4Y','  C4Z',' C4X3',' C4Y3',' C4Z3',
     &             '  C2a','  C2b','  C2c','  C2d','  C2e','  C2f',
     &             '  C3Z',' C3Z2','  C6Z',' C6Z5','  C2A','  C2B',
     &             '  C2C','  C2D','    I',' IC2X',' IC2Y',' IC2Z',
     &             ' IC3l',' IC3m',' IC3n',' IC3o','IC3l2','IC3m2',
     &             'IC3n2','IC3o2',' IC4X',' IC4Y',' IC4Z','IC4X3',
     &             'IC4Y3','IC4Z3',' IC2a',' IC2b',' IC2c',' IC2d',
     &             ' IC2e',' IC2f',' IC3Z','IC3Z2',' IC6Z','IC6Z5',
     &             ' IC2A',' IC2B',' IC2C',' IC2D'/
C
      DATA ANG_ROT /'PI/4','PI/2','PI/6','PI/3'/
C
      DATA AT_ADD /'        ','<--- NEW'/
C
      DATA SMALL /0.001/
C
      PRINT 444
C
      DO JS=1,64
        IS_WRONG(JS)=0
      ENDDO
      IROT=0
      NAT_NEW=NATCLU
C
      IF(I_GR.EQ.2) THEN
        I_SB=1
      ELSE
        I_SB=0
      ENDIF
C
 111  DO JAT=1,NAT_NEW
        NEQAT(JAT)=0
      ENDDO
C
C    Calculates the distance between the atoms and all
C    their neighbors
C
      DIST0(1)=0.
      DO JAT1=1,NAT_NEW
        INEW_AT(JAT1)=0
        NABOR=0
        DO JAT2=1,NAT_NEW
          IF(JAT1.EQ.JAT2) GO TO 10
          NABOR=NABOR+1
          NBRZ(NABOR,JAT1)=NZAT(JAT2)
          RAB=SQRT((COORD(1,JAT1)-COORD(1,JAT2))*(COORD(1,JAT1)-
     &    COORD(1,JAT2))+(COORD(2,JAT1)-COORD(2,JAT2))*(COORD(2,JAT1)-
     &    COORD(2,JAT2))+(COORD(3,JAT1)-COORD(3,JAT2))*(COORD(3,JAT1)-
     &    COORD(3,JAT2)))
          IF((JAT2.GT.JAT1).AND.(RAB.LT.SMALL)) GOTO 895
          IF(JAT1.EQ.1) DIST0(JAT2)=RAB
          DNBR(NABOR,JAT1)=RAB
 10       CONTINUE
        ENDDO
      ENDDO
C
C   Generates the set of class distances to the absorber
C
      CALL ORDRE(NAT_NEW,DIST0,NDIST,DIST1)
C
C     Determines the prototypical and equivalent atoms:
C     two atoms are equivalent if they have the same type,
C     they are at the same distance from the emitter
C     and they have the same geometrical environment.
C     This part of the routine has been adapted from
C     R. Gunnella and C. R. Natoli.
C
      JTYP_C_M=0
      NEQAT_M=0
      JTYP=2
      NASHUFF=2
      NABORS=NAT_NEW-1
C
      NATYP(1)=1
      NQAT(1)=1
      NCORR(1,1)=1
      NEQAT(1)=1
      NCHTYP(1)=1
      I_Z(1,1)=1
      Z_L(1,1)=1.
      Z_M1(1,1)=ONEC
      Z_M2(1,1)=ONEC
      SYM_AT(1,1)=COORD(1,1)
      SYM_AT(2,1)=COORD(2,1)
      SYM_AT(3,1)=COORD(3,1)
      NZ_AT(1)=NZAT(1)
      CHEM(1)=CHEM_OLD(1)
      INEW_AT(1)=0
      N_NEW_OLD(1)=1
C
      NATYP(2)=1
      NQAT(2)=2
      NCORR(1,2)=2
      NEQAT(2)=2
      NCHTYP(2)=2
      I_Z(1,2)=1
      Z_L(1,2)=1.
      Z_M1(1,2)=ONEC
      Z_M2(1,2)=ONEC
      SYM_AT(1,2)=COORD(1,2)
      SYM_AT(2,2)=COORD(2,2)
      SYM_AT(3,2)=COORD(3,2)
      NZ_AT(2)=NZAT(2)
      CHEM(2)=CHEM_OLD(2)
      INEW_AT(2)=0
      N_NEW_OLD(2)=2
C
      JNUM=0
      DO JAT1=2,NAT_NEW
        IF(JNUM.GT.NEQAT_M) THEN
          NEQAT_M=JNUM
          JTYP_C_M=JTYP
        ENDIF
        JNUM=1
        NA1P1=JAT1+1
        IF(NA1P1.GT.NAT_NEW) GO TO 32
        DO JAT2=NA1P1,NAT_NEW
          IF(NZAT(JAT1).NE.NZAT(JAT2)) GO TO 30
          N1=JAT1-1
          N2=JAT2-1
          IF(ABS(DNBR(N1,1)-DNBR(N2,1)).GT.SMALL) GOTO 30
          DO NAB=1,NABORS
            IUSED(NAB)=0
          ENDDO
          DO NABOR1=1,NABORS
            NZT=NBRZ(NABOR1,JAT1)
            RABT=DNBR(NABOR1,JAT1)
            EQUIV=.FALSE.
            DO NABOR2=1,NABORS
              IF(IUSED(NABOR2).EQ.1) GOTO 22
              IF(NBRZ(NABOR2,JAT2).NE.NZT) GOTO 22
              IF(ABS(DNBR(NABOR2,JAT2)-RABT).GT.SMALL) GOTO 22
              EQUIV=.TRUE.
              IUSED(NABOR2)=1
              GOTO 23
 22           CONTINUE
            ENDDO
 23         IF(.NOT.EQUIV) GO TO 30
          ENDDO
          IF(NEQAT(JAT1).EQ.0) THEN
            JTYP=JTYP+1
            NASHUFF=NASHUFF+1
            NCORR(JNUM,JTYP)=NASHUFF
            I_Z(JNUM,JTYP)=1
            Z_L(JNUM,JTYP)=1.
            Z_M1(JNUM,JTYP)=ONEC
            Z_M2(JNUM,JTYP)=ONEC
            SYM_AT(1,NASHUFF)=COORD(1,JAT1)
            SYM_AT(2,NASHUFF)=COORD(2,JAT1)
            SYM_AT(3,NASHUFF)=COORD(3,JAT1)
            N_NEW_OLD(NASHUFF)=JAT1
            NEQAT(JAT1)=JTYP
            NCHTYP(JTYP)=NTYP(JAT1)
            NZ_AT(NASHUFF)=NZAT(JAT1)
            CHEM(NASHUFF)=CHEM_OLD(JAT1)
            NQAT(NASHUFF)=JTYP
            IF(JAT1.GT.NATCLU) INEW_AT(NASHUFF)=1
          ENDIF
          IF(NEQAT(JAT2).EQ.0) THEN
            JNUM=JNUM+1
            NATYP(JTYP)=JNUM
            NASHUFF=NASHUFF+1
            NCORR(JNUM,JTYP)=NASHUFF
            SYM_AT(1,NASHUFF)=COORD(1,JAT2)
            SYM_AT(2,NASHUFF)=COORD(2,JAT2)
            SYM_AT(3,NASHUFF)=COORD(3,JAT2)
            N_NEW_OLD(NASHUFF)=JAT2
            NEQAT(JAT2)=JTYP
            NCHTYP(JTYP)=NTYP(JAT2)
            NZ_AT(NASHUFF)=NZAT(JAT2)
            CHEM(NASHUFF)=CHEM_OLD(JAT2)
            NQAT(NASHUFF)=JTYP
            IF(JAT2.GT.NATCLU) INEW_AT(NASHUFF)=1
         ENDIF
30       CONTINUE
        ENDDO
32      IF(NEQAT(JAT1).EQ.0) THEN
          JTYP=JTYP+1
          NATYP(JTYP)=JNUM
          NASHUFF=NASHUFF+1
          NCORR(JNUM,JTYP)=NASHUFF
          I_Z(JNUM,JTYP)=1
          Z_L(JNUM,JTYP)=1.
          Z_M1(JNUM,JTYP)=ONEC
          Z_M2(JNUM,JTYP)=ONEC
          SYM_AT(1,NASHUFF)=COORD(1,JAT1)
          SYM_AT(2,NASHUFF)=COORD(2,JAT1)
          SYM_AT(3,NASHUFF)=COORD(3,JAT1)
          N_NEW_OLD(NASHUFF)=JAT1
          NEQAT(JAT1)=JTYP
          NCHTYP(JTYP)=NTYP(JAT1)
          NZ_AT(NASHUFF)=NZAT(JAT1)
          CHEM(NASHUFF)=CHEM_OLD(JAT1)
          NQAT(NASHUFF)=JTYP
          IF(JAT1.GT.NATCLU) INEW_AT(NASHUFF)=1
        ENDIF
        CONTINUE
      ENDDO
      N_PROT=JTYP
C
C     Stop if the maximal number of prototypical and equivalent
C                atoms are not correctly dimensionned
C
      IF(N_PROT.GT.NATP_M) THEN
        WRITE(IUO1,897) N_PROT
        STOP
      ENDIF
      IF(NEQAT_M.GT.NAT_EQ_M) THEN
        WRITE(IUO1,898) NEQAT_M
        STOP
      ENDIF
C
      DO JAT=1,NAT_NEW
        COORD1(1,JAT)=SYM_AT(1,JAT)-SYM_AT(1,1)
        COORD1(2,JAT)=SYM_AT(2,JAT)-SYM_AT(2,1)
        COORD1(3,JAT)=SYM_AT(3,JAT)-SYM_AT(3,1)
      ENDDO
C
C     Test of all symmetry operations for the largest
C              symmetry class
C
      NROTCLUS=0
 556  NSYM_M(1)=1
      NSIZE_C=1
      DO JROT=2,64
        DO JNUM1=1,NEQAT_M
          JAT1=NCORR(JNUM1,JTYP_C_M)
          MATCH=.FALSE.
          SYM_AT1(1,JAT1)=S_M(1,1,JROT)*COORD1(1,JAT1)+S_M(1,2,
     &    JROT)*COORD1(2,JAT1)+S_M(1,3,JROT)*COORD1(3,JAT1)
          SYM_AT1(2,JAT1)=S_M(2,1,JROT)*COORD1(1,JAT1)+S_M(2,2,
     &    JROT)*COORD1(2,JAT1)+S_M(2,3,JROT)*COORD1(3,JAT1)
          SYM_AT1(3,JAT1)=S_M(3,1,JROT)*COORD1(1,JAT1)+S_M(3,2,
     &    JROT)*COORD1(2,JAT1)+S_M(3,3,JROT)*COORD1(3,JAT1)
          DO JNUM2=1,NEQAT_M
            JAT2=NCORR(JNUM2,JTYP_C_M)
            AD=ABS(COORD1(1,JAT2)-SYM_AT1(1,JAT1))+ABS(COORD1(2,
     &      JAT2)-SYM_AT1(2,JAT1))+ABS(COORD1(3,JAT2)-SYM_AT1(3,JAT1))
            IF(AD.LT.SMALL) THEN
              MATCH=.TRUE.
            ENDIF
          ENDDO
          IF(.NOT.MATCH) GOTO 333
        ENDDO
        NSIZE_C=NSIZE_C+1
        NSYM_M(NSIZE_C)=JROT
 333    CONTINUE
      ENDDO
C
C      Test on all classes of the symmetry operations that work
C            for the largest class
C
      NSYM_G(1)=1
      NSIZE_GR=1
      I_THREE=0
      I_INCRG=0
      I_INVG=0
      DO JSYM=2,NSIZE_C
        JROT=NSYM_M(JSYM)
        DO JTYP=1,N_PROT
          NEQATS=NATYP(JTYP)
          DO JNUM1=1,NEQATS
            JAT1=NCORR(JNUM1,JTYP)
            MATCH=.FALSE.
            SYM_AT1(1,JAT1)=S_M(1,1,JROT)*COORD1(1,JAT1)+S_M(1,2,
     &      JROT)*COORD1(2,JAT1)+S_M(1,3,JROT)*COORD1(3,JAT1)
            SYM_AT1(2,JAT1)=S_M(2,1,JROT)*COORD1(1,JAT1)+S_M(2,2,
     &      JROT)*COORD1(2,JAT1)+S_M(2,3,JROT)*COORD1(3,JAT1)
            SYM_AT1(3,JAT1)=S_M(3,1,JROT)*COORD1(1,JAT1)+S_M(3,2,
     &      JROT)*COORD1(2,JAT1)+S_M(3,3,JROT)*COORD1(3,JAT1)
            DO JNUM2=1,NEQATS
              JAT2=NCORR(JNUM2,JTYP)
              AD=ABS(COORD1(1,JAT2)-SYM_AT1(1,JAT1))+ABS(
     &        COORD1(2,JAT2)-SYM_AT1(2,JAT1))+ABS(COORD1(3,JAT2)-
     &        SYM_AT1(3,JAT1))
              IF(AD.LT.SMALL) THEN
                MATCH=.TRUE.
              ENDIF
            ENDDO
            IF(.NOT.MATCH) GOTO 335
          ENDDO
        ENDDO
        IF(MATCH) THEN
          NSIZE_GR=NSIZE_GR+1
          NSYM_G(NSIZE_GR)=JROT
          IF(JROT.EQ.25) I_THREE=1
          IF(JROT.EQ.33) I_INVG=1
        ENDIF
 335    CONTINUE
      ENDDO
      IF(NSIZE_GR.GT.48) GOTO 998
C
C      Set up of the parameter used to check if a larger
C         group containing the inversion I can be constructed
C         if the inversion is physically acceptable as an approximation
C
      IF((I_GR.EQ.2).AND.(I_INVG.EQ.0)) THEN
        IF(I_INV.EQ.1) THEN
          I_INCRG=1
        ENDIF
      ENDIF
C
      IF(NROTCLUS.LE.2) THEN
        IF(NSIZE_GR.GT.1) THEN
          WRITE(IUO1,699) NSIZE_GR
        ELSE
          WRITE(IUO1,698) NSIZE_GR
        ENDIF
        IF(NSIZE_GR.EQ.1) THEN
          WRITE(IUO1,705) NAME_S(NSYM_G(1))
        ELSEIF(NSIZE_GR.EQ.2) THEN
          WRITE(IUO1,704) NAME_S(NSYM_G(1)),NAME_S(NSYM_G(2))
        ELSEIF(NSIZE_GR.EQ.3) THEN
          WRITE(IUO1,703) NAME_S(NSYM_G(1)),NAME_S(NSYM_G(2)),
     &    NAME_S(NSYM_G(3))
        ELSEIF(NSIZE_GR.EQ.4) THEN
          WRITE(IUO1,702) NAME_S(NSYM_G(1)),NAME_S(NSYM_G(2)),
     &    NAME_S(NSYM_G(3)),NAME_S(NSYM_G(4))
        ELSEIF(NSIZE_GR.EQ.6) THEN
          WRITE(IUO1,701) (NAME_S(NSYM_G(JROT)), JROT=1,6)
        ELSEIF(NSIZE_GR.GE.8) THEN
          WRITE(IUO1,700) (NAME_S(NSYM_G(JROT)), JROT=1,NSIZE_GR)
        ENDIF
        IF(NROTCLUS.GT.0) THEN
          IF(I_THREE.EQ.0) THEN
            WRITE(IUO1,706) ANG_ROT(NROTCLUS)
          ELSE
            WRITE(IUO1,706) ANG_ROT(NROTCLUS+2)
          ENDIF
        ENDIF
      ENDIF
C
C      Finding the cluster's symmetry group
C
 555  JGROUP=0
      IF(NSIZE_GR.EQ.1) THEN
        JGROUP=1
      ELSEIF(NSIZE_GR.EQ.2) THEN
        IF(NSYM_G(2).EQ.33) THEN
          JGROUP=2
        ELSEIF(NSYM_G(2).EQ.4) THEN
          JGROUP=3
        ELSEIF(NSYM_G(2).EQ.36) THEN
          JGROUP=4
        ENDIF
      ELSEIF(NSIZE_GR.EQ.3) THEN
        JGROUP=16
      ELSEIF(NSIZE_GR.EQ.4) THEN
        IF(NSYM_G(3).EQ.33) THEN
          JGROUP=5
        ELSEIF((NSYM_G(3).EQ.3).OR.(NSYM_G(3).EQ.19)) THEN
          JGROUP=6
        ELSEIF((NSYM_G(3).EQ.34).OR.(NSYM_G(3).EQ.51)) THEN
          JGROUP=7
        ELSEIF(NSYM_G(3).EQ.15) THEN
          JGROUP=9
        ELSEIF(NSYM_G(3).EQ.47) THEN
          JGROUP=10
        ENDIF
      ELSEIF(NSIZE_GR.EQ.6) THEN
        IF(NSYM_G(4).EQ.26) THEN
          JGROUP=18
        ELSEIF((NSYM_G(4).EQ.34).OR.(NSYM_G(4).EQ.35)) THEN
          JGROUP=19
        ELSEIF(NSYM_G(4).EQ.33) THEN
          JGROUP=17
        ELSEIF(NSYM_G(4).EQ.26) THEN
          JGROUP=21
        ELSEIF(NSYM_G(4).EQ.36) THEN
          JGROUP=22
        ENDIF
      ELSEIF(NSIZE_GR.EQ.8) THEN
        IF(NSYM_G(4).EQ.33) THEN
          IF(NSYM_G(8).EQ.50) THEN
            JGROUP=11
          ELSE
            JGROUP=8
          ENDIF
        ELSE
          IF(NSYM_G(5).EQ.15) THEN
            JGROUP=12
          ELSE
            IF(NSYM_G(3).EQ.15) THEN
              JGROUP=13
            ELSEIF((NSYM_G(3).EQ.3).OR.(NSYM_G(3).EQ.19)) THEN
              JGROUP=14
            ENDIF
          ENDIF
        ENDIF
      ELSEIF(NSIZE_GR.EQ.12) THEN
        IF(NSYM_G(5).EQ.5) THEN
          JGROUP=28
        ELSEIF(NSYM_G(7).EQ.33) THEN
          IF(NSYM_G(12).EQ.60) THEN
            JGROUP=23
          ELSE
            JGROUP=20
          ENDIF
        ELSE
          IF(NSYM_G(3).EQ.3) THEN
            JGROUP=24
          ELSEIF(NSYM_G(9).EQ.59) THEN
            JGROUP=26
          ELSEIF(NSYM_G(5).EQ.27) THEN
            JGROUP=25
          ENDIF
        ENDIF
      ELSEIF(NSIZE_GR.EQ.16) THEN
        JGROUP=15
      ELSEIF(NSIZE_GR.EQ.24) THEN
        IF(NSYM_G(17).EQ.57) THEN
          JGROUP=27
        ELSEIF(NSYM_G(17).EQ.17) THEN
          JGROUP=31
        ELSEIF(NSYM_G(17).EQ.49) THEN
          JGROUP=30
        ELSEIF(NSYM_G(17).EQ.37) THEN
          JGROUP=29
        ENDIF
      ELSEIF(NSIZE_GR.EQ.48) THEN
        JGROUP=32
      ENDIF
c      IF((JGROUP.GT.0).AND.(NEQAT_M.LE.NSIZE_GR)) THEN
      IF(JGROUP.GT.0) THEN
        WRITE(IUO1,886)
        WRITE(IUO1,880) NAME_G(JGROUP)
      ELSE
C
C  If no group is found, the cluster is rotated by pi/4 ,or pi/6 for
C    hexagonal structures (i. e. when C3z=25 is present), around Oz
C    at most twice to account for possible misnaming of the Ox axis.
C    Then, the search for the point-group is restarted.
C
        IF(NROTCLUS.LT.2) THEN
          NROTCLUS=NROTCLUS+1
          DO JPROT=1,N_PROT
            NEQATS=NATYP(JPROT)
            DO JNUM=1,NEQATS
              JAT=NCORR(JNUM,JPROT)
              X=COORD1(1,JAT)
              Y=COORD1(2,JAT)
              IF(I_THREE.EQ.0) THEN
                COORD1(1,JAT)=0.707107*(X-Y)
                COORD1(2,JAT)=0.707107*(X+Y)
              ELSE
                COORD1(1,JAT)=0.866025*X-0.500000*Y
                COORD1(2,JAT)=0.500000*X+0.866025*Y
              ENDIF
              COORD1(3,JAT)=COORD1(3,JAT)
              SYM_AT(1,JAT)=COORD1(1,JAT)+SYM_AT(1,1)
              SYM_AT(2,JAT)=COORD1(2,JAT)+SYM_AT(2,1)
              SYM_AT(3,JAT)=COORD1(3,JAT)+SYM_AT(3,1)
            ENDDO
          ENDDO
          GOTO 556
        ELSE
          IF((JGROUP.EQ.0).AND.(NROTCLUS.EQ.2))  WRITE(IUO1,881)
        ENDIF
      ENDIF
      IF(JGROUP.GE.28) WRITE(IUO1,697)
      IF((JGROUP.GT.0) .AND.(I_INCRG.EQ.1)) WRITE(IUO1,722)
C
C  Recovery of the original cluster when no group has been found
C
      IF((NROTCLUS.EQ.2).AND.(JGROUP.EQ.0)) THEN
        DO JPROT=1,N_PROT
          NEQATS=NATYP(JPROT)
          DO JNUM=1,NEQATS
            JAT=NCORR(JNUM,JPROT)
            X=COORD1(1,JAT)
            Y=COORD1(2,JAT)
            IF(I_THREE.EQ.0) THEN
              COORD1(1,JAT)=Y
              COORD1(2,JAT)=-X
            ELSE
              COORD1(1,JAT)=0.500000*X+0.866025*Y
              COORD1(2,JAT)=-0.866025*X+0.500000*Y
            ENDIF
            COORD1(3,JAT)=COORD1(3,JAT)
            SYM_AT(1,JAT)=COORD1(1,JAT)+SYM_AT(1,1)
            SYM_AT(2,JAT)=COORD1(2,JAT)+SYM_AT(2,1)
            SYM_AT(3,JAT)=COORD1(3,JAT)+SYM_AT(3,1)
          ENDDO
        ENDDO
        NROTCLUS=NROTCLUS+1
        GOTO 556
      ENDIF
C
C  If still no group is found, or if the group can be augmented by
C   the inversion (as an approximation), check of the other symmetries that
C   should be present to account for group properties (i. e. those
C   obtained by multiplication of those that are effectively present)
C
      IF((JGROUP.EQ.0).OR.(I_INCRG.EQ.1)) THEN
        IF(I_INCRG.EQ.1) THEN
          NSIZE_GR=NSIZE_GR+1
          NSYM_G(NSIZE_GR)=33
        ENDIF
        NSIZE_RE=NSIZE_GR
 553    I_NEW=0
        DO JSYM1=2,NSIZE_RE
          JS1=NSYM_G(JSYM1)
          DO JSYM2=JSYM1,NSIZE_RE
            JS2=NSYM_G(JSYM2)
            DO I=1,3
              DO J=1,3
                S_MUL(I,J)=0.
                DO K=1,3
                  S_MUL(I,J)=S_MUL(I,J)+S_M(I,K,JS1)*S_M(K,
     &            J,JS2)
                ENDDO
              ENDDO
            ENDDO
            I_EQ=0
            DO JSYM3=1,NSIZE_RE
              JS3=NSYM_G(JSYM3)
              I_OLD=0
              DO I=1,3
                DO J=1,3
                  S1=S_MUL(I,J)
                  S2=S_M(I,J,JS3)
                  IF(ABS(S1-S2).LT.SMALL) I_OLD=I_OLD+1
                ENDDO
              ENDDO
              IF(I_OLD.EQ.9) I_EQ=1
            ENDDO
            IF(I_EQ.EQ.0) THEN
              I_NEW=I_NEW+1
              J_RE=NSIZE_GR+I_NEW
              DO JS4=2,64
                I_OLD=0
                DO I=1,3
                  DO J=1,3
                    S1=S_MUL(I,J)
                    S2=S_M(I,J,JS4)
                    IF(ABS(S1-S2).LT.SMALL) I_OLD=I_OLD+1
                  ENDDO
                ENDDO
                IF(I_OLD.EQ.9) THEN
                  NSYM_G(J_RE)=JS4
                  NSIZE_RE=NSIZE_RE+1
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        IF(I_NEW.GT.0) GOTO 553
C
        IF(NSIZE_RE.GT.NSIZE_GR) THEN
          WRITE(IUO1,696) NSIZE_RE
          IF(NSIZE_RE.EQ.2) THEN
            WRITE(IUO1,704) NAME_S(NSYM_G(1)),NAME_S(NSYM_G(2))
          ELSEIF(NSIZE_RE.EQ.3) THEN
            WRITE(IUO1,703) NAME_S(NSYM_G(1)),NAME_S(NSYM_G(2)),
     &      NAME_S(NSYM_G(3))
          ELSEIF(NSIZE_RE.EQ.4) THEN
            WRITE(IUO1,702) NAME_S(NSYM_G(1)),NAME_S(NSYM_G(2)),
     &      NAME_S(NSYM_G(3)),NAME_S(NSYM_G(4))
          ELSEIF(NSIZE_RE.EQ.6) THEN
            WRITE(IUO1,701) (NAME_S(NSYM_G(JROT)), JROT=1,6)
          ELSEIF(NSIZE_RE.GE.8) THEN
              WRITE(IUO1,700) (NAME_S(NSYM_G(JROT)), JROT=1,
     &        NSIZE_RE)
          ENDIF
        ENDIF
      ENDIF
C
C  If no group has been found, it means that the cluster does not
C   retain the whole symmetry of the crystal. When I_GR = 1 or 2,
C   the largest group consistent with the symmetry operations
C   found is sought. Then a new cluster is built to support
C             all the symmetries of this group.
C
      I_BACK=0
      IF(((JGROUP.EQ.0).OR.(I_INCRG.EQ.1)).AND.(I_GR.GE.1)) THEN
C
C  Search for the different point-groups containing the NSIZE_RE
C   symmetries found. If the cluster can not support the
C   inversion I (I_INV=0), a test is made on the content
C   of the point-groups to suppress those that contain I.
C
        WRITE(IUO1,707)
        WRITE(IUO1,708)
C
C  Input cluster
C
        NB_GR=0
        DO JG=1,32
          IF((I_INV.EQ.0).AND.(I_33(JG).EQ.1)) GOTO 223
          JS_MAX=SIZE_G(JG)
          IF(JS_MAX.LT.NSIZE_RE) GOTO 223
C
          ICHECK_S=0
          DO JSYM=1,NSIZE_RE
            JS_RE=NSYM_G(JSYM)
C
            DO JS_GR1=1,JS_MAX
              JS_NE1=CONT_G(JS_GR1,JG)
              IWRONG=0
              IF(IROT.GT.0) THEN
                DO JS=1,IROT
                  IF(JS_NE1.EQ.IS_WRONG(JS)) THEN
                    IWRONG=1
                  ENDIF
                ENDDO
              ENDIF
              IF(JS_NE1.EQ.JS_RE) THEN
                ICHECK_S=ICHECK_S+1
              ENDIF
            ENDDO
C
          ENDDO
C
          IF(IWRONG.EQ.1) GOTO 223
          IF(ICHECK_S.EQ.NSIZE_RE) THEN
            NB_GR=NB_GR+1
            OLD_G(NB_GR)=JG
          ENDIF
 223      CONTINUE
C
        ENDDO
        NB_GR1=NB_GR
        WRITE(IUO1,709) (NAME_G(OLD_G(J)), J=1,NB_GR1)
        WRITE(IUO1,710)
C
C  Input cluster rotated
C
        DO JG=1,32
          IF((I_INV.EQ.0).AND.(I_33(JG).EQ.1)) GOTO 225
          IOLD=0
          DO J_OLD=1,NB_GR1
            IF(OLD_G(J_OLD).EQ.JG) THEN
               IOLD=IOLD+1
            ENDIF
          ENDDO
          IF(IOLD.NE.0) GOTO 225
          JS_MAX=SIZE_G(JG)
          IF(JS_MAX.LT.NSIZE_RE) GOTO 225
C
          ICHECK_S=0
          DO JSYM=1,NSIZE_RE
            JS_RE=NSYM_G(JSYM)
C
            DO JS_GR1=1,JS_MAX
              JS_GR2=49-JS_GR1
              JS_NE2=CONT_G(JS_GR2,JG)
              IWRONG=0
              IF(IROT.GT.0) THEN
                DO JS=1,IROT
                  IF(JS_NE2.EQ.IS_WRONG(JS)) THEN
                    IWRONG=1
                  ENDIF
                ENDDO
              ENDIF
              IF(JS_NE2.EQ.JS_RE) THEN
                ICHECK_S=ICHECK_S+1
              ENDIF
            ENDDO
C
          ENDDO
C
          IF(IWRONG.EQ.1) GOTO 225
          IF(ICHECK_S.EQ.NSIZE_RE) THEN
            NB_GR=NB_GR+1
            OLD_G(NB_GR)=JG
          ENDIF
 225      CONTINUE
C
        ENDDO
        IF(NB_GR1.LT.NB_GR) THEN
          WRITE(IUO1,713)
          WRITE(IUO1,711) (NAME_G(OLD_G(J)), J=NB_GR1+1,NB_GR)
          WRITE(IUO1,712)
        ENDIF
        IF(I_INV.EQ.0) THEN
          IF(I_GR.EQ.2) THEN
            WRITE(IUO1,721)
          ELSE
            WRITE(IUO1,726)
          ENDIF
        ENDIF
        CALL ORDRE2(NB_GR,OLD_G,NB_G,NEW_G)
C
C  Construction of the augmented cluster consistent with the group JGROUP_M.
C    Note that this group must be consistent with the original cluster.
C    This is checked through 3 criteria :
C
C         1. No atom must be generated above the surface plane
C
C         2. If an atom generated by a symmetry of JGROUP_M coincides
C             with one of the original cluster, it must have the same
C             atomic number
C
C         3. Every new atom generated by a symmetry of JGROUP_M must
C             be outside the original cluster except those meeting
C             the previous criterion
C
C  When one of this criteria is not satisfied, it means that the group
C   JGROUP_M can not accomodate the original cluster. Hence, a smaller
C   group must be looked for
C
C  An extra criterion is used when I_GR = 1 :
C
C         4. No surface atom can be transformed into a bulk atom
C              and likewise no bulk atom into a surface atom
C
C
        JG_MIN=NEW_G(1)
        I_MINUS=0
 222    JGROUP_M=NEW_G(NB_G-I_MINUS)
        IF(JGROUP_M.LT.JG_MIN) THEN
          IF(I_GR.GE.1) THEN
            I_GR=0
            WRITE(IUO1,723)
            GOTO 111
          ELSE
            WRITE(IUO1,719)
            STOP
          ENDIF
        ENDIF
        JS_M=SIZE_G(JGROUP_M)
        WRITE(IUO1,714) NAME_G(JGROUP_M)
        I_END=0
        IROT=0
        NAT=NAT_NEW
C
        DO JS=2,JS_M
          JROT=CONT_G(JS,JGROUP_M)
C
          DO J_AT=2,NAT
            IF(I_END.EQ.1) GOTO 224
            X_NEW=S_M(1,1,JROT)*(COORD(1,J_AT)-COORD(1,1)) + S_M(
     &      1,2,JROT)*(COORD(2,J_AT)-COORD(2,1)) +S_M(1,3,JROT)*(COORD(
     &      3,J_AT)-COORD(3,1)) + COORD(1,1)
            Y_NEW=S_M(2,1,JROT)*(COORD(1,J_AT)-COORD(1,1)) + S_M(
     &      2,2,JROT)*(COORD(2,J_AT)-COORD(2,1)) +S_M(2,3,JROT)*(COORD(
     &      3,J_AT)-COORD(3,1)) + COORD(2,1)
            Z_NEW=S_M(3,1,JROT)*(COORD(1,J_AT)-COORD(1,1)) + S_M(
     &      3,2,JROT)*(COORD(2,J_AT)-COORD(2,1)) +S_M(3,3,JROT)*(COORD(
     &      3,J_AT)-COORD(3,1)) + COORD(3,1)
C
C  Check for criterion 1
C
            IF(Z_NEW.GT.VALZ_MAX) THEN
              WRITE(IUO1,715) (COORD(J,J_AT), J=1,3),X_NEW,
     &        Y_NEW,Z_NEW
              WRITE(IUO1,716) NAME_S(JROT),VALZ_MAX
              IROT=IROT+1
              IS_WRONG(IROT)=JROT
            ENDIF
            IF(IROT.GT.0) THEN
              I_END=1
              GOTO 224
            ENDIF
            NZ_NEW=NZAT(J_AT)
            I_OLD=0
            I_IN_OK=0
C
C  Check for criterion 2
C
            DO J_AT_P=2,NAT
              D2=(X_NEW-COORD(1,J_AT_P))*(X_NEW-COORD(1,J_AT_P)
     &        ) + (Y_NEW-COORD(2,J_AT_P))*(Y_NEW-COORD(2,J_AT_P)) +(
     &        Z_NEW-COORD(3,J_AT_P))*(Z_NEW-COORD(3,J_AT_P))
              NZ_OLD=NZAT(J_AT_P)
              IF(D2.LT.SMALL) THEN
                I_OLD=I_OLD+1
                IF(NZ_NEW.NE.NZ_OLD) THEN
                  IROT=IROT+1
                  IS_WRONG(IROT)=JROT
                  WRITE(IUO1,715) (COORD(J,J_AT), J=1,3),
     &            X_NEW,Y_NEW,Z_NEW
                  WRITE(IUO1,717)  CHEM_OLD(J_AT),J_AT,
     &            NAME_S(JROT),CHEM_OLD(J_AT_P),J_AT_P
                ELSE
                  I_IN_OK=1
                ENDIF
              ENDIF
            ENDDO
            IF(I_IN_OK.EQ.1) GOTO 226
C
C  Check for criterion 3
C
            DO JPLAN=1,NPLAN
              IF(ABS(Z_NEW-Z_PLAN(JPLAN)).LT.SMALL) THEN
                IF(X_NEW.GT.X_MIN(JPLAN)) THEN
                  IF(X_NEW.LT.X_MAX(JPLAN)) THEN
                    IF(Y_NEW.GT.Y_MIN(JPLAN)) THEN
                      IF(Y_NEW.LT.Y_MAX(JPLAN)) THEN
                        IROT=IROT+1
                        IS_WRONG(IROT)=JROT
                        WRITE(IUO1,715) (COORD(J,
     &                  J_AT), J=1,3),X_NEW,Y_NEW,Z_NEW
                        WRITE(IUO1,720) NAME_S(JROT)
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
C
C  Check for criterion 4
C
           IF(I_SB.EQ.0) THEN
            Z_DIF=ABS(Z_NEW-COORD(3,J_AT))
            Z_TYP=ABS(Z_NEW-VALZ_MAX)
            IF(Z_DIF.GT.SMALL) THEN
              WRITE(IUO1,715) (COORD(J,J_AT), J=1,3),X_NEW,
     &        Y_NEW,Z_NEW
              IF(Z_TYP.LT.SMALL) THEN
                WRITE(IUO1,725) NAME_S(JROT),I_GR
              ELSE
                WRITE(IUO1,724) NAME_S(JROT),I_GR
              ENDIF
              IROT=IROT+1
              IS_WRONG(IROT)=JROT
            ENDIF
            IF(IROT.GT.0) THEN
              I_END=1
              GOTO 224
            ENDIF
           ENDIF
C
 226       IF(I_OLD.EQ.0) THEN
              NAT=NAT+1
              IF(NAT.GT.NATCLU_M) THEN
                WRITE(IUO1,718) NAT
                STOP
              ENDIF
              COORD(1,NAT)=X_NEW
              COORD(2,NAT)=Y_NEW
              COORD(3,NAT)=Z_NEW
              VALZ(NAT)=Z_NEW
              CHEM_OLD(NAT)=CHEM_OLD(J_AT)
              NZAT(NAT)=NZAT(J_AT)
              NTYP(NAT)=NTYP(J_AT)
            ENDIF
 224        CONTINUE
          ENDDO
C
        ENDDO
C
        I_BACK=1
        IF(IROT.GT.0) THEN
          I_MINUS=I_MINUS+1
          GOTO 222
        ENDIF
      ENDIF
      IF(I_BACK.EQ.1) THEN
        NAT_NEW=NAT
        GOTO 111
      ENDIF
C
C  Writes the classes of atoms by increasing distance
C
      WRITE(IUO1,888)
      DO JDIST=NDIST,1,-1
        DO JTYP=1,N_PROT
          NMAX=NATYP(JTYP)
          DO JCLASS=1,NMAX
            N=NCORR(JCLASS,JTYP)
            D1=SQRT(COORD1(1,N)*COORD1(1,N)+COORD1(2,N)*COORD1(2,
     &      N)+COORD1(3,N)*COORD1(3,N))
            IF(ABS(D1-DIST1(JDIST)).LT.SMALL) THEN
              WRITE(IUO1,557) N,SYM_AT(1,N),SYM_AT(2,N),SYM_AT(
     &        3,N),NQAT(N),JCLASS,CHEM(N),DIST1(JDIST),AT_ADD(INEW_AT(
     &        N)+1)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
C  Writes the augmented symmetrized cluster into the output
C          file OUTFILE4 for further use if necessary
C
      OPEN(UNIT=IUO4, FILE=OUTFILE4, STATUS='UNKNOWN')
      WRITE(IUO4,778) IPHA
      DO JTYP=1,N_PROT
        NMAX=NATYP(JTYP)
        DO JCLASS=1,NMAX
          N=NCORR(JCLASS,JTYP)
          X=SYM_AT(1,N)/CUNIT
          Y=SYM_AT(2,N)/CUNIT
          Z=SYM_AT(3,N)/CUNIT
          IF(IPHA.EQ.0) THEN
            WRITE(IUO4,125) N,CHEM(N),NZ_AT(N),X,Y,Z,JTYP
          ELSEIF(IPHA.EQ.1) THEN
            WRITE(IUO4,779) N,CHEM(N),NZ_AT(N),X,Y,Z,JTYP
          ELSEIF(IPHA.EQ.2) THEN
            WRITE(IUO4,*) NZ_AT(N),X,Y,Z,JTYP
          ENDIF
        ENDDO
      ENDDO
      CLOSE(IUO4)
C
C  Construction of the symmetry operations list.
C    Associates the appropriate symmetry operation
C         for all equivalent atoms
C
      IF(IPRINT.GE.2) WRITE(IUO1,887)
 438  CONTINUE
      I_CALC_ROT=0
      DO JTYP=1,N_PROT
        ISYM(1,JTYP)=1
        NEQATS=NATYP(JTYP)
        IF(NEQATS.EQ.1) THEN
          IF(JTYP.GT.N_PROT) GOTO 555
          GOTO 338
        ENDIF
        NAP=NCORR(1,JTYP)
        DO JNUM=1,NEQATS
	  IF(JNUM.GE.2) THEN
	    JSYM_0=2
	  ELSE
	    JSYM_0=1
	  ENDIF
          NA=NCORR(JNUM,JTYP)
          X=COORD1(1,NAP)
          Y=COORD1(2,NAP)
          Z=COORD1(3,NAP)
          DO JSYM=JSYM_0,NSIZE_GR
            JROT=NSYM_G(JSYM)
            SYM_AT1(1,NAP)=S_M(1,1,JROT)*X+S_M(1,2,JROT)*Y+S_M(1,
     &      3,JROT)*Z
            SYM_AT1(2,NAP)=S_M(2,1,JROT)*X+S_M(2,2,JROT)*Y+S_M(2,
     &      3,JROT)*Z
            SYM_AT1(3,NAP)=S_M(3,1,JROT)*X+S_M(3,2,JROT)*Y+S_M(3,
     &      3,JROT)*Z
            AD=ABS(COORD1(1,NA)-SYM_AT1(1,NAP))+ABS(COORD1(2,NA)-
     &      SYM_AT1(2,NAP))+ABS(COORD1(3,NA)-SYM_AT1(3,NAP))
            IF(AD.LT.SMALL)THEN
              ISYM(JNUM,JTYP)=JROT
              I_Z(JNUM,JTYP)=IZ(JROT)
              Z_L(JNUM,JTYP)=ZL(JROT)
              Z_M1(JNUM,JTYP)=ZM1(JROT)
              Z_M2(JNUM,JTYP)=ZM2(JROT)
              IF(IZ(JROT).EQ.0) THEN
                I_CALC_ROT=I_CALC_ROT+1
              ENDIF
              GOTO 404
            ENDIF
          ENDDO
          IF(ISYM(JNUM,JTYP).EQ.0) THEN
            N_PROT=N_PROT+1
            NCHTYP(N_PROT)=NCHTYP(JTYP)
            NCORR(1,N_PROT)=NA
            NATYP(N_PROT)=1
            I_Z(1,N_PROT)=1
            Z_L(1,N_PROT)=1.
            Z_M1(1,N_PROT)=ONEC
            Z_M2(1,N_PROT)=ONEC
            DO JCHANGE=JNUM,NEQATS-1
              NCORR(JCHANGE,JTYP)=NCORR(JCHANGE+1,JTYP)
            ENDDO
            NATYP(JTYP)=NATYP(JTYP)-1
            GOTO 438
          ENDIF
 404      CONTINUE
          IF((IPRINT.GE.2).AND.(NSIZE_GR.GT.1)) THEN
	    JR=ISYM(JNUM,JTYP)
            WRITE(IUO1,849) JTYP,JNUM,NCORR(JNUM,JTYP),NAME_S(JR)
     &      ,JR,Z_L(JNUM,JTYP),Z_M1(JNUM,JTYP),Z_M2(JNUM,JTYP),I_Z(
     &      JNUM,JTYP)
          ENDIF
        ENDDO
	WRITE(IUO1,*) '     '
 338    CONTINUE
      ENDDO
C
      GAIN_G=REAL(NAT_NEW)/REAL(N_PROT)
      WRITE(IUO1,854) GAIN_G
C
C     Test of the symmetry operations leaving the prototypical
C       atoms invariant. Associates the apropriate symmetry
C          relation and/or selection rule for each atom.
C
C            NAT_SYM(J) is the number of prototypical atoms in the
C                       various symmetry sets :
C
C                             J = 1 : atom 0
C                             J = 2 : z axis
C                             J = 3 : x0y plane
C                             J = 4 : other atoms
C
      IF(IPRINT.GE.2) WRITE(IUO1,889)
      NAT_SYM(1)=1
      NAT_SYM(2)=0
      NAT_SYM(3)=0
      NAT_SYM(4)=0
      I_SET(1)=1
C
C   Loop on the prototypical atoms
C
      DO JTYP=1,N_PROT
C
        ISTEP_L(JTYP)=1
        ISTEP_M(JTYP)=1
        I_REL_MP(JTYP)=0
        I_LM(JTYP)=0
        I_Z_P(JTYP)=0
        Z_L_P(JTYP)=1.
        Z_M_P(JTYP)=ONEC
C
        NSYM_P=1
        NSYM_PT=1
        INV_P(JTYP,NSYM_P)=1
        INV_PT(JTYP,NSYM_PT)=1
C
        JAT_P=NCORR(1,JTYP)
        X=COORD1(1,JAT_P)
        Y=COORD1(2,JAT_P)
        Z=COORD1(3,JAT_P)
        X_A=ABS(X)
        Y_A=ABS(Y)
        Z_A=ABS(Z)
        IF(JTYP.GT.1) THEN
          IF((X_A+Y_A.LT.SMALL).AND.(Z_A.GE.SMALL)) THEN
            NAT_SYM(2)=NAT_SYM(2)+1
            I_SET(JTYP)=2
          ELSEIF((Z_A.LT.SMALL).AND.(X_A+Y_A.GE.SMALL)) THEN
            NAT_SYM(3)=NAT_SYM(3)+1
            I_SET(JTYP)=3
          ELSEIF(((X_A+Y_A).GE.SMALL).AND.(Z_A.GE.SMALL)) THEN
            NAT_SYM(4)=NAT_SYM(4)+1
            I_SET(JTYP)=4
          ENDIF
        ENDIF
C
C     Loop on the symmetries keeping the cluster unchanged
C
        DO JSYM=2,NSIZE_GR
C
          JROT=NSYM_G(JSYM)
          X1=S_M(1,1,JROT)*X+S_M(1,2,JROT)*Y+S_M(1,3,JROT)*Z
          Y1=S_M(2,1,JROT)*X+S_M(2,2,JROT)*Y+S_M(2,3,JROT)*Z
          Z1=S_M(3,1,JROT)*X+S_M(3,2,JROT)*Y+S_M(3,3,JROT)*Z
          AD=ABS(X-X1)+ABS(Y-Y1)+ABS(Z-Z1)
C
C       Case of an atom invariant by the symmetry JROT
C
          IF(AD.LT.SMALL) THEN
            NSYM_PT=NSYM_PT+1
            INV_PT(JTYP,NSYM_PT)=JROT
            IF(IZ(JROT).NE.0) THEN
              I_Z_P(JTYP)=IZ(JROT)
              NSYM_P=NSYM_P+1
              INV_P(JTYP,NSYM_P)=JROT
              ISL=ISTEP_L(JTYP)
              ISM=ISTEP_M(JTYP)
C
C       Case of an atom off the z axis
C
              IF((ABS(X).GE.SMALL).OR.(ABS(Y).GE.SMALL)) THEN
C
C         Symmetry = IC2z
C
                IF(JROT.EQ.36) THEN
                  ISTEP_M(JTYP)=MAX(ISM,2)
                  I_LM(JTYP)=1
C
C         Symmetry = C2u or IC2u_perp
C
                ELSE
                  I_REL_MP(JTYP)=1
                  Z_L_P(JTYP)=ZL(JROT)
                  Z_M_P(JTYP)=ZM1(JROT)
                ENDIF
C
C       Case of an atom on the z axis but different from the absorber
C
              ELSE
                IF(ABS(Z).GE.SMALL) THEN
C
C         Symmetry = C2z
C
                  IF(JROT.EQ.4) THEN
                    ISTEP_M(JTYP)=MAX(ISM,2)
C
C         Symmetry = C4z
C
                  ELSEIF(JROT.EQ.15) THEN
                    ISTEP_M(JTYP)=MAX(ISM,4)
C
C         Symmetry = C4z3
C
                  ELSEIF(JROT.EQ.18) THEN
                    ISTEP_M(JTYP)=MAX(ISM,4)
C
C         Symmetry = C3z
C
                  ELSEIF(JROT.EQ.25) THEN
                    ISTEP_M(JTYP)=MAX(ISM,3)
C
C         Symmetry = C3z2
C
                  ELSEIF(JROT.EQ.26) THEN
                    ISTEP_M(JTYP)=MAX(ISM,3)
C
C         Symmetry = C6z
C
                  ELSEIF(JROT.EQ.27) THEN
                    ISTEP_M(JTYP)=MAX(ISM,6)
C
C         Symmetry = C6z5
C
                  ELSEIF(JROT.EQ.28) THEN
                    ISTEP_M(JTYP)=MAX(ISM,6)
C
C         Symmetry = IC2u
C
                  ELSEIF(JROT.GT.33) THEN
                    I_REL_MP(JTYP)=1
                    I_Z_P(JTYP)=IZ(JROT)
                    Z_L_P(JTYP)=ZL(JROT)
                    Z_M_P(JTYP)=ZM1(JROT)
                  ENDIF
C
C       Case of  atom 0 (the absorber)
C
                ELSE
C
C         Symmetry = C2z or IC2z
C
                  IF((JROT.EQ.4).OR.(JROT.EQ.36)) THEN
                    ISTEP_M(JTYP)=MAX(ISM,2)
                    IF(JROT.EQ.36) THEN
                      I_LM(JTYP)=1
                    ENDIF
C
C         Symmetry = C4z or IC4z
C
                  ELSEIF((JROT.EQ.15).OR.(JROT.EQ.47)) THEN
                    ISTEP_M(JTYP)=MAX(ISM,4)
                    IF(JROT.EQ.47) THEN
                      I_LM(JTYP)=1
                    ENDIF
C
C         Symmetry = C4z3 or IC4z3
C
                  ELSEIF((JROT.EQ.18).OR.(JROT.EQ.50)) THEN
                    ISTEP_M(JTYP)=MAX(ISM,4)
                    IF(JROT.EQ.50) THEN
                      I_LM(JTYP)=1
                    ENDIF
C
C         Symmetry = C3z or IC3z
C
                  ELSEIF((JROT.EQ.25).OR.(JROT.EQ.57)) THEN
                    ISTEP_M(JTYP)=MAX(ISM,3)
                    IF(JROT.EQ.57) THEN
                      ISTEP_L(JTYP)=MAX(ISL,2)
                    ENDIF
C
C         Symmetry = C3z2 or IC3z2
C
                  ELSEIF((JROT.EQ.26).OR.(JROT.EQ.58)) THEN
                    ISTEP_M(JTYP)=MAX(ISM,3)
                    IF(JROT.EQ.58) THEN
                      ISTEP_L(JTYP)=MAX(ISL,2)
                    ENDIF
C
C         Symmetry = C6z or IC6z
C
                  ELSEIF((JROT.EQ.27).OR.(JROT.EQ.59)) THEN
                    ISTEP_M(JTYP)=MAX(ISM,6)
                    IF(JROT.EQ.59) THEN
                      I_LM(JTYP)=1
                    ENDIF
C
C         Symmetry = C6z5 or IC6z5
C
                  ELSEIF((JROT.EQ.28).OR.(JROT.EQ.60)) THEN
                    ISTEP_M(JTYP)=MAX(ISM,6)
                    IF(JROT.EQ.60) THEN
                      I_LM(JTYP)=1
                    ENDIF
C
C         Symmetry = I
C
                  ELSEIF(JROT.EQ.33) THEN
                    ISTEP_L(JTYP)=MAX(ISL,2)
C
C         Symmetry = C2u or IC2u_perp with u within (x0y)
C
                  ELSE
                    IF(IZ(JROT).EQ.-1) THEN
                      I_REL_MP(JTYP)=1
                      Z_L_P(JTYP)=ZL(JROT)
                      Z_M_P(JTYP)=ZM1(JROT)
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
C
C   Finding the symmetry group (if any) associated to each prototypical atom
C
        JGROUP=0
        IF(NSYM_PT.EQ.1) THEN
          JGROUP=1
        ELSEIF(NSYM_PT.EQ.2) THEN
          IF(INV_PT(JTYP,2).EQ.33) THEN
            JGROUP=2
          ELSEIF(INV_PT(JTYP,2).EQ.4) THEN
            JGROUP=3
          ELSEIF(INV_PT(JTYP,2).EQ.36) THEN
            JGROUP=4
          ENDIF
        ELSEIF(NSYM_PT.EQ.3) THEN
          JGROUP=16
        ELSEIF(NSYM_PT.EQ.4) THEN
          IF(INV_PT(JTYP,3).EQ.33) THEN
            JGROUP=5
          ELSEIF((INV_PT(JTYP,3).EQ.3).OR.(INV_PT(JTYP,3).EQ.19)) 
     &    THEN
            JGROUP=6
          ELSEIF((INV_PT(JTYP,3).EQ.34).OR.(INV_PT(JTYP,3).EQ.51)) 
     &    THEN
            JGROUP=7
          ELSEIF(INV_PT(JTYP,3).EQ.15) THEN
            JGROUP=9
          ELSEIF(INV_PT(JTYP,3).EQ.47) THEN
            JGROUP=10
          ENDIF
        ELSEIF(NSYM_PT.EQ.6) THEN
          IF(INV_PT(JTYP,4).EQ.26) THEN
            JGROUP=18
          ELSEIF((INV_PT(JTYP,4).EQ.34).OR.(INV_PT(JTYP,4).EQ.35)) 
     &    THEN
            JGROUP=19
          ELSEIF(INV_PT(JTYP,4).EQ.33) THEN
            JGROUP=17
          ELSEIF(INV_PT(JTYP,4).EQ.26) THEN
            JGROUP=21
          ELSEIF(INV_PT(JTYP,4).EQ.36) THEN
            JGROUP=22
          ENDIF
        ELSEIF(NSYM_PT.EQ.8) THEN
          IF(INV_PT(JTYP,4).EQ.33) THEN
            IF(INV_PT(JTYP,8).EQ.50) THEN
              JGROUP=11
            ELSE
              JGROUP=8
            ENDIF
          ELSE
            IF(INV_PT(JTYP,5).EQ.15) THEN
              JGROUP=12
            ELSE
              IF(INV_PT(JTYP,3).EQ.15) THEN
                JGROUP=13
              ELSEIF((INV_PT(JTYP,3).EQ.3).OR.(INV_PT(JTYP,3).
     &        EQ.19)) THEN
                JGROUP=14
              ENDIF
            ENDIF
           ENDIF
        ELSEIF(NSYM_PT.EQ.12) THEN
          IF(INV_PT(JTYP,5).EQ.5) THEN
            JGROUP=28
          ELSEIF(INV_PT(JTYP,7).EQ.33) THEN
            IF(INV_PT(JTYP,12).EQ.60) THEN
              JGROUP=23
            ELSE
              JGROUP=20
            ENDIF
          ELSE
            IF(INV_PT(JTYP,3).EQ.3) THEN
              JGROUP=24
            ELSEIF(INV_PT(JTYP,9).EQ.59) THEN
              JGROUP=26
            ELSEIF(INV_PT(JTYP,5).EQ.27) THEN
              JGROUP=25
            ENDIF
          ENDIF
        ELSEIF(NSYM_PT.EQ.16) THEN
          JGROUP=15
        ELSEIF(NSYM_PT.EQ.24) THEN
          IF(INV_PT(JTYP,17).EQ.57) THEN
            JGROUP=27
          ELSEIF(INV_PT(JTYP,17).EQ.17) THEN
            JGROUP=31
          ELSEIF(INV_PT(JTYP,17).EQ.49) THEN
            JGROUP=30
          ELSEIF(INV_PT(JTYP,17).EQ.37) THEN
            JGROUP=29
          ENDIF
        ELSEIF(NSYM_PT.EQ.48) THEN
          JGROUP=32
        ENDIF
        GR(JTYP)=JGROUP
C
        IF((IPRINT.GE.2).AND.(NSIZE_GR.GT.1)) THEN
          WRITE(IUO1,851) JTYP,JAT_P,SYM_AT(1,JAT_P),SYM_AT(2,
     &    JAT_P),SYM_AT(3,JAT_P),I_SET(JTYP),NSYM_PT,NSYM_P,NAME_G(GR(
     &    JTYP)),(NAME_S(INV_P(JTYP,JS)),JS=1,NSYM_P)
          ENDIF
      ENDDO
      WRITE(IUO1,852)
      GAIN_B=0.
      DO JTYP=1,N_PROT
        NGAIN_B(JTYP)=ISTEP_L(JTYP)*ISTEP_M(JTYP)*(I_REL_MP(JTYP)+1)
        GAIN_B=GAIN_B+NGAIN_B(JTYP)
        WRITE(IUO1,853) JTYP,I_Z_P(JTYP),INT(Z_L_P(JTYP)),Z_M_P(JTYP)
     &  ,ISTEP_L(JTYP),ISTEP_M(JTYP),I_LM(JTYP),I_REL_MP(JTYP),NGAIN_B(
     &  JTYP)
      ENDDO
      GAIN_B=GAIN_B/FLOAT(N_PROT)
      WRITE(IUO1,855) GAIN_B
C
C  Calculation and storage of r^{l}_{m m'}(pi/2) for the specific
C            cubic group symmetry operations
C
      IF(I_CALC_ROT.GT.0) THEN
        PIS2=1.570796
        CALL DJMN(PIS2,R_PIS2,LI_M+1)
      ENDIF
C
C  Construction of the inverse matrices used for point-groups
C
      DO I=1,3
        DO J=1,3
          DO JOP=1,NSIZE_GR
            JSYM=NSYM_G(JOP)
            S_INV(JSYM,J,I)=S_M(I,J,JSYM)
          ENDDO
        ENDDO
      ENDDO
      GOTO 999
C
 895  WRITE(IUO1,896) JAT1,JAT2
      STOP
 998  WRITE(IUO1,997) NSIZE_GR
      STOP
 999  PRINT 445
C
C  Input/output formats
C
 125  FORMAT(2X,I4,5X,A2,5X,I2,3F10.4,12X,I4)
 444  FORMAT(////,5X,'++++++++++++++++++++  SYMMETRIZATION OF THE ',
     &'CLUSTER   +++++++++++++++++++',/)
 445  FORMAT(//,5X,'+++++++++++++++++++++++++++++++++++++++++++++++','+
     &++++++++++++++++++++++++',//)
 557  FORMAT(10X,I3,3X,'(',F7.3,',',F7.3,',',F7.3,')',3X,I3,3X,I3,3X,
     &A2,3X,F7.3,2X,A8)
 696  FORMAT(///,6X,I2,' SYMMETRY OPERATIONS AT LEAST SHOULD BE ',
     &'PRESENT IN THE CLUSTER :',/)
 697  FORMAT(///,16X,' THE THREE-FOLD AXES HAVE BEEN CODED AS : ',//,
     &30X,'alpha  ----->  l',/,30X,'beta   ----->  m',/,30X,'gamma  ---
     &-->  n',/,30X,'delta  ----->  o')
 698  FORMAT(///,17X,I1,' SYMMETRY OPERATION FOUND IN THE CLUSTER :',/)
 699  FORMAT(///,16X,I2,' SYMMETRY OPERATIONS FOUND IN THE CLUSTER :',/
     &)
 700  FORMAT(12X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,/,12X,A5,
     &2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,/,     12X,A5,2X,A5,2X,
     &A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,/,     12X,A5,2X,A5,2X,A5,2X,A5,
     &2X,A5,2X,A5,2X,A5,2X,A5,/,     12X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,
     &A5,2X,A5,2X,A5,/,     12X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,
     &2X,A5,/,     12X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,/, 
     &    12X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5)
 701  FORMAT(19X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5)
 702  FORMAT(26X,A5,2X,A5,2X,A5,2X,A5)
 703  FORMAT(30X,A5,2X,A5,2X,A5)
 704  FORMAT(33X,A5,2X,A5)
 705  FORMAT(37X,A5)
 706  FORMAT(/,19X,'AFTER ROTATION OF THIS CLUSTER BY ',A4,/)
 707  FORMAT(//,'    ----------',9X,'CONSTRUCTION OF A NEW CLUSTER :',
     &8X,'----------')
 708  FORMAT(//,10X,'THE DIFFERENT GROUPS THAT COULD SUPPORT ','THE 
     &CLUSTER ARE :',/)
 709  FORMAT(28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,3X),
     &/,28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,3X),/,
     &28X,4(A3,3X))
 710  FORMAT(/,28X,'FOR THE INPUT CLUSTER ')
 711  FORMAT(/,28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,
     &3X),/,28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,3X),/
     &,28X,4(A3,3X))
 712  FORMAT(/,24X,'FOR THE ROTATED INPUT CLUSTER')
 713  FORMAT(//,36X,'AND :')
 714  FORMAT(//,25X,'---> TRYING THE GROUP : ',A3)
 715  FORMAT(/,13X,'GROUP IMPOSSIBLE : ATOM (',F7.3,',',F7.3,',',F7.3,
     &')',/,19X,' TRANSFORMED INTO (',F7.3,',',F7.3,',',F7.3,')')
 716  FORMAT(20X,'THE NEW ATOM GENERATED BY ',A5,' IS ABOVE',/,20X,
     &'THE SURFACE PLANE LOCATED AT Z = ',F7.3)
 717  FORMAT(27X,'NEW ATOM OF TYPE ',A2,' (No ',I4,') GENERATED BY ',
     &A5,/,27X,'AT THE POSITION OF AN ATOM OF TYPE ',A2,' (No ',I4,')',
     &/,27X,'IN THE ORIGINAL CLUSTER')
 718  FORMAT(///,5X,'<<<<<<<<<<  NATCLU_M TOO SMALL IN THE .inc FILE ',
     &'TO INCREASE  >>>>>>>>>>',/,5X,'<<<<<<<<<<        THE ','CLUSTER.
     & SHOULD BE AT LEAST ',I4,'       >>>>>>>>>>')
 719  FORMAT(///,3X,'<<<<<<<<<<  ERROR : NO GROUP WAS FOUND TO ',
     &'ACCOMODATE THE CLUSTER  >>>>>>>>>>')
 720  FORMAT(20X,'THE NEW ATOM GENERATED BY ',A5,' IS INSIDE',/,25X,
     &'THE ORIGINAL CLUSTER, BUT NOT',/,29X,'AT AN ATOMIC POSITION')
 721  FORMAT(//,11X,' THE INVERSION I IS NOT CONSISTENT WITH THIS ',
     &'CLUSTER',/,17X,'AS THE ABSORBER IS CLOSER TO THE SURFACE',/,21X,
     &'THAN TO THE BOTTOM OF THE CLUSTER')
 722  FORMAT(//,16X,'THE CLUSTER CAN ACCOMODATE THE INVERSION I',//,
     &13X,'---> BUILDING A LARGER CLUSTER INVARIANT BY I')
 723  FORMAT(//,4X,'--------   IMPOSSIBLE TO AUGMENT THE CLUSTER TO ',
     &'SUPPORT I    --------',//,15X,'---> RESTARTING WITH THE 
     &ORIGINAL CLUSTER ...')
 724  FORMAT(20X,'THE NEW ATOM GENERATED BY ',A5,' IS OF BULK TYPE',/,
     &25X,'WHILE THE ORIGINAL ONE IS OF SURFACE TYPE',/,29X,'---> 
     &IMPOSSIBLE WITH I_GR = ',I1)
 725  FORMAT(20X,'THE NEW ATOM GENERATED BY ',A5,' IS OF SURFACE TYPE',
     &/,25X,'WHILE THE ORIGINAL ONE IS OF BULK TYPE',/,29X,'---> 
     &IMPOSSIBLE WITH I_GR = ',I1)
 726  FORMAT(//,11X,' THE INVERSION I IS NOT CONSISTENT WITH THIS ',
     &'CALCULATION',/,10X,'AS SURFACE AND BULK ATOMS HAVE TO ','BE 
     &DISCRIMINATED (I_GR=1)')
 778  FORMAT(30X,I1)
 779  FORMAT(2X,I4,5X,A2,5X,I2,3F10.4,I5)
 849  FORMAT(1X,I4,4X,I2,3X,I4,2X,A5,' = ',I2,2X,F6.3,2X,'(',F6.3,',',
     &F6.3,')',2X,'(',F6.3,',',F6.3,')',4X,I2)
 851  FORMAT(1X,I4,3X,I4,2X,'(',F7.3,',',F7.3,',',F7.3,')',2X,I1,3X,I2,
     &1X,I2,2X,A3,2X,4(1X,A5),/,57X,4(1X,A5),/,57X,4(1X,A5),/,57X,4(1X,
     &A5),/,57X,4(1X,A5),/,57X,4(1X,A5),/,57X,4(1X,A5),/,57X,4(1X,A5),/
     &,57X,4(1X,A5),/,57X,4(1X,A5),/,57X,4(1X,A5),/,57X,4(1X,A5),/,57X)
 852  FORMAT(///,10X,'SELECTION RULE AND/OR RELATION ON THE MATRIX',' 
     &ELEMENTS OF TAU :',//,4X,'CLASS',2X,'I_Z',2X,'Z_L',8X,'Z_M',9X,
     &'ISTEP_L',2X,'ISTEP_M',2X,'I_LM',2X,'I_REL_MP',2X,'GAIN',/)
 853  FORMAT(4X,I3,4X,I2,3X,I2,3X,'(',F6.3,',',F6.3,')',6X,I1,8X,I1,7X,
     &I1,7X,I1,6X,I2)
 854  FORMAT(//,16X,'-----> EXPECTED GAIN FOR THE GLOBAL LEVEL : ',F5.
     &2,//)
 855  FORMAT(///,18X,'-----> EXPECTED GAIN FOR THE BASIS LEVEL : ',F5.
     &2,//)
 880  FORMAT(33X,A3)
 881  FORMAT(/,'    ----------  THESE OPERATIONS DON''T FORM A',' 
     &SYMMETRY GROUP  ----------', /'    ----------  i. e. THE CLUSTER 
     &IS NOT CORRECTLY',' TRUNCATED  ----------',/,'    ----------     
     &  FROM THE POINT OF VIEW OF ','SYMMETRY       ----------')
 886  FORMAT(///,25X,'CLUSTER SYMMETRY GROUP :',/,33X,A4)
 887  FORMAT(///,10X,'SYMMETRY OPERATIONS ASSOCIATED WITH THE ',
     &'EQUIVALENT ATOMS : ',//,2X,'CLASS',2X,'ATOM',2X,' No ',3X,
     &'SYMMETRY',3X,'Z**L',8X,'Z**M',13X,'Z**M''',7X,'DELTA',/)
 888  FORMAT(///,20X,'CONTENTS OF THE SYMMETRIZED CLUSTER',/,18X,'BY 
     &INCREASING DISTANCE TO THE ABSORBER :',//,11X,'No',13X,'(X,Y,Z)',
     &11X,'CLASS',1X,'ATOM',1X,'ChSp',4X,'DIST',//)
 889  FORMAT(///,10X,'SYMMETRY OPERATIONS LEAVING THE PROTOTYPICAL ',
     &'ATOMS INVARIANT : ',//,27X,'G = 1  ----->  atom 0 only',/,27X,
     &'G = 2  ----->  atom along 0z',/,27X,'G = 3  ----->  atom within 
     &x0y',/,27X,'G = 4  ----->  other atom',//,19X,'ST : total number 
     &of symmetries leaving',/,29X,'the prototypical atom invariant',//
     &,19X,'SE : number of symmetries taken into account',/,29X,'(
     &Euler angle BETA = 0 or pi)',//,2X,'CLASS',2X,' No ',11X,'(X,Y,Z)
     &',10X,'G',3X,'ST',1X,'SE',1X,'GROUP',7X,'SE SYMMETRIES',/)
 896  FORMAT(///,'<<<<<<<<<<  ERROR IN THE COORDINATES OF THE  ATOMS',
     &'  >>>>>>>>>>',/,'<<<<<<<<<<  ATOMS ',I4,' AND ',I4,' ARE 
     &IDENTICAL  >>>>>>>>>>')
 897  FORMAT(///,11X,'<<<<<<<<<<  NATP_M TOO SMALL IN THE INCLUDE ',
     &'FILE  >>>>>>>>>>',/,11X,'<<<<<<<<<<        SHOULD BE',' AT 
     &LEAST ',I4,'         >>>>>>>>>>')
 898  FORMAT(///,12X,'<<<<<<<<<<  NAT_EQ_M TOO SMALL IN THE INCLUDE',' 
     &FILE  >>>>>>>>>>',/,12X,'<<<<<<<<<< SHOULD BE AT LEAST ',I4,'  
     &>>>>>>>>>>')
 997  FORMAT(///,'<<<<<<<<<<  ',I2,' SYMMETRIES HAVE BEEN FOUND. THIS 
     &','>>>>>>>>>>',/,'<<<<<<<<<<  EXCEEDS THE SIZE OF THE ','LARGEST 
     &POINT-GROUP  >>>>>>>>>>')
C
      RETURN
C
      END
