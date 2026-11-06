memalloc_src               := memalloc/dim_mod.f memalloc/modules.f memalloc/allocation.f 
cluster_gen_src            := $(sort $(wildcard cluster_gen/*.f))
common_sub_src             := $(sort $(wildcard common_sub/*.f))
renormalization_src        := $(sort $(wildcard renormalization/*.f))
eig_common_src             := $(sort $(wildcard eig/common/*.f))
eig_pw_src                 := $(sort $(wildcard eig/pw/*.f))

SRCS   = $(memalloc_src) $(cluster_gen_src) $(common_sub_src) $(renormalization_src) $(eig_common_src) $(eig_pw_src)
MAIN_F = eig/pw/main.f
SO     = _eig_pw.so

include ../../../options.mk
