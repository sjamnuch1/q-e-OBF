# Makefile for simple_ip

include ../../make.inc

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../Modules

OBFDIAGOBJS = \
obf_diag_objects.o \
start_end.o \
interpolation_obf.o \
input_obf_diag.o \
tetra_mod1.o \
diagonalization.o \
read_igwx.o \
read_wfc.o \
write_wfc.o \
unk_from_wfc.o \
dielectric.o  \
dielectric_spin.o


QEMODS = $(BASEMODS)
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a

LIBMIN= 

TLDEPS=pwlibs

all : tldeps obf_diag.x  libobf_diag.a

obf_diag.x : obf_diag.o $(OBFDIAGOBJS) $(PWOBJS) $(QEMODS)  $(LIBMIN) 
	$(MPIF90) $(LDFLAGS) -o $@ \
		obf_diag.o $(OBFDIAGOBJS) $(PWOBJS) $(QEMODS) $(LIBMIN) $(QELIBS) 
	- ( cd ../../bin; ln -fs ../OBF/obf_diag/$@ . )

libobf_diag.a : $(OBFDIAGOBJS)
	$(AR) $(ARFLAGS) $@ $(OBFDIAGOBJS)

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f -v simple_ip.x  *.o *~ *.F90 *.d *.mod *.i work.pc

include make.depend

# DO NOT DELETE
