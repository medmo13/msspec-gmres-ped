memalloc_src               := memalloc/dim_mod.f memalloc/modules.f memalloc/allocation.f 
cluster_gen_src            := $(sort $(wildcard cluster_gen/*.f))
common_sub_src             := $(sort $(wildcard common_sub/*.f))
renormalization_src        := $(sort $(wildcard renormalization/*.f))
phd_se_noso_nosp_nosym_src := $(sort $(wildcard phd_se_noso_nosp_nosym/*.f))

SRCS   = $(memalloc_src) $(cluster_gen_src) $(common_sub_src) $(renormalization_src) $(phd_se_noso_nosp_nosym_src)
MAIN_F = phd_se_noso_nosp_nosym/main.f
SO     = _phd_se_noso_nosp_nosym.so

include ../../../options.mk
