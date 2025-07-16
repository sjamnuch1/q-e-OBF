# Makefile for  obf

include ../../make.inc
#LIBMIN = ../minpack/minpacklib.a

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../obf_ham \
	  $(MOD_FLAG)../../Modules/

OBFOBJS = \
  stop_pp.o \
  init_us_2_max.o \
  read_export.o \
  kpoint_mesh.o \
  openfile_ham.o \
  read_kpoint.o \
  obf_ham.o \
  input_ham.o \
  khamiltonian.o \
  khamiltonian_byfile.o \
  khamiltonian_spin.o \
  khamiltonian_byfile_spin.o \

QEMODS = $(BASEMODS)
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a 
MODOBJ = ../../Modules/libqemod.a

TLDEPS= pwlibs

all : tldeps obf_ham.x

obf_ham.x : obf_ham.o libham.a  $(OBFOBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
		obf_ham.o libham.a $(PWOBJS) $(QEMODS) $(QELIBS) $(LIBMIN)
	- ( cd ../../bin ; ln -fs ../OBF/obf_ham/$@ . )

tldeps :
	if test -n "$(TLDEPS)" ; then   ( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi


libham.a : $(OBFOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@


clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L libobf.a

include make.depend
# DO NOT DELETE

