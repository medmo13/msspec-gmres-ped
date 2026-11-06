memalloc_src               := memalloc/dim_mod.f memalloc/modules.f memalloc/allocation.f 
cluster_gen_src            := $(sort $(wildcard cluster_gen/*.f))
common_sub_src             := $(sort $(wildcard common_sub/*.f))
renormalization_src        := $(sort $(wildcard renormalization/*.f))
#eig_common_src             := $(sort $(wildcard eig/common/*.f))
eig_common_src             := $(filter-out eig/common/lapack_eig.f, $(sort $(wildcard eig/common/*.f)))
eig_mi_src                 := $(sort $(wildcard eig/mi/*.f))

SRCS   = $(memalloc_src) $(cluster_gen_src) $(common_sub_src) $(renormalization_src) $(eig_common_src) $(eig_mi_src)
MAIN_F = eig/mi/main.f
SO     = _eig_mi.so

include ../../../options.mk
