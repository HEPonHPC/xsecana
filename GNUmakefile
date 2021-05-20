PACKAGE     := XSecAna
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
SUBDIRS     := CAFAna

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_eigen.mk


LIBLINK := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)

# Explicitly link in root libraries. Why is this necessary? Shouldn't
# arch_spec_root.mk do that? Nonetheless, without this we get undefined symbol
# errors...
override LIBLIBS += `root-config --libs` -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -L$(IFDHC_FQ_DIR)/lib/ -lifdh

