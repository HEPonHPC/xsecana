
BINS =  test_compile \
	test_hist \
	test_interface \
	test_systematic \
	test_simple_efficiency \
	test_simple_flux \
	test_simple_signal_estimator \
	test_simple_xsec 	
SIMPLEBINS=$(BINS)

LIB_TYPE    := 

include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_cafanacore.mk
include SoftRelTools/arch_spec_srproxy.mk
include SoftRelTools/arch_spec_eigen.mk
include SoftRelTools/arch_spec_stan.mk

override BINLIBS += `root-config --libs` \
	-L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) \
	-L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) \
	-L$(IFDHC_FQ_DIR)/lib/ \
	-L$(CAFANACORE_LIB) \
	-lifdh \
	-lStandardRecord \
	-lStandardRecord_dict \
	-lStandardRecordProxy \
	-lCAFAnaCoreExt \
	-lCAFAnaCore \
	-lXSecAnaCAFAna

