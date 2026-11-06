


      MODULE DIM_MOD
        IMPLICIT NONE
        INTEGER NATP_M, NATCLU_M, NAT_EQ_M, N_CL_L_M
        INTEGER NE_M, NL_M
        INTEGER LI_M, NEMET_M
        INTEGER NO_ST_M
        INTEGER NDIF_M
        INTEGER NSO_M
        INTEGER NTEMP_M
        INTEGER NODES_EX_M
C
        INTEGER NSPIN_M
C
        INTEGER NTH_M, NPH_M
        INTEGER NDIM_M
        INTEGER N_TILT_M
        INTEGER N_ORD_M
C
        INTEGER NGR_M
        INTEGER NPATH_M

C       ===============================================================
        INTEGER NLP_M, NLA_M
        INTEGER N_MU_M, N_NU_M
        INTEGER NATM
        INTEGER LINMAX, LINMAXA
        INTEGER LINFMAX
        INTEGER NLAMBDA_M
        INTEGER NSPIN2_M
        INTEGER NT_M
        INTEGER NCG_M
        INTEGER N_BESS, N_GAUNT
        INTEGER NLTWO
        INTEGER NLMM
C       ===============================================================
      CONTAINS
        SUBROUTINE INIT_DIM()
          NLP_M=NL_M
          NLA_M=NL_M
        
          N_MU_M=NO_ST_M
          N_NU_M=NO_ST_M/2
    
          NATM=NATP_M+3
          
          LINMAX=NLP_M*NLP_M
          LINMAXA=NLA_M*NLA_M
        
          LINFMAX=(LI_M+2)*(LI_M+2)
        
          NLAMBDA_M=(NO_ST_M+2)*(NO_ST_M+1)/2
        
          NSPIN2_M=3*NSPIN_M-2
        
          NT_M=(NL_M-1)*(1+(NSPIN_M-1)*NL_M)
        
          NCG_M=4*LI_M+2

C          N_BESS=100*NL_M
C          N_GAUNT=5*NL_M
          N_BESS=300*NL_M
          N_GAUNT=10*NL_M

          NLTWO=2*NL_M
          NLMM=LINMAX*NGR_M
        END SUBROUTINE INIT_DIM
      END MODULE DIM_MOD
