.PHONY: all phd_se phd_mi phd_gm phd_ce eig_mi eig_pw comp_curve clean

all: phd_se phd_mi phd_gm phd_ce eig_mi eig_pw comp_curve

phd_se:
	@+$(MAKE) -f phd_se_noso_nosp_nosym.mk all

phd_mi:
	@+$(MAKE) -f phd_mi_noso_nosp_nosym.mk all

phd_gm:
	@+$(MAKE) -f phd_gm_noso_nosp_nosym.mk all

phd_ce:
	@+$(MAKE) -f phd_ce_noso_nosp_nosym.mk all

eig_mi:
	@+$(MAKE) -f eig_mi.mk all

eig_pw:
	@+$(MAKE) -f eig_pw.mk all

comp_curve:
	@+$(MAKE) -f comp_curve.mk all

clean::
	@+$(MAKE) -f phd_se_noso_nosp_nosym.mk $@
	@+$(MAKE) -f phd_mi_noso_nosp_nosym.mk $@
	@+$(MAKE) -f phd_gm_noso_nosp_nosym.mk $@
	@+$(MAKE) -f phd_ce_noso_nosp_nosym.mk $@
	@+$(MAKE) -f eig_mi.mk $@
	@+$(MAKE) -f eig_pw.mk $@
	@+$(MAKE) -f comp_curve.mk $@
