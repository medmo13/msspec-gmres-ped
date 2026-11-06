memalloc_src               := memalloc/dim_mod.f memalloc/modules.f memalloc/allocation.f 
cluster_gen_src            := $(sort $(wildcard cluster_gen/*.f))
common_sub_src             := $(sort $(wildcard common_sub/*.f))
renormalization_src        := $(sort $(wildcard renormalization/*.f))
phd_gm_noso_nosp_nosym_src := $(filter-out phd_gm_noso_nosp_nosym/lapack_axb.f, $(sort $(wildcard phd_gm_noso_nosp_nosym/*.f)))

SRCS   = $(memalloc_src) $(cluster_gen_src) $(common_sub_src) $(renormalization_src) $(phd_gm_noso_nosp_nosym_src)
MAIN_F = phd_gm_noso_nosp_nosym/main.f
SO     = _phd_gm_noso_nosp_nosym.so

include ../../../options.mk
F2PYFLAGS += -lgomp

# -ffast-math implies -ffinite-math-only, which removes the IEEE_IS_FINITE guards of the solver (D44): keep them for this file.
$(BUILDDIR)/phd_gm_noso_nosp_nosym/axb_mat_ms_la.o: FFLAGS += -fno-finite-math-only
