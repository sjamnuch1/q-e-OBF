include ../make.inc

default: all

all: obf_basisa

obf_basisa :
	( cd obf_basis ; $(MAKE) all || exit 1 )

obf_hama :
	( cd obf_ham ; $(MAKE) all || exit 1 )

obf_diaga :
	( cd obf_diag ; $(MAKE) all || exit 1 )
