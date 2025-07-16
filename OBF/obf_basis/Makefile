# Makefile for  obf

include ../../make.inc
#LIBMIN = ../minpack/minpacklib.a

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../obf_basis \
	  $(MOD_FLAG)../../Modules/

OBFOBJS = \
  stop_pp.o \
  read_export.o \
  read_kpoint.o \
  data_structure_obf.o \
  openfile_basis.o \
  input_basis.o \
  wfc_basis.o \
  s_basis.o \
  f_basis.o \
  read_obf.o \
  split_obf.o \
  gk_l2gmap.o \
  gk_sort_simple.o \
  gk_sort_count.o \
  gk_sort_ngw_lmax.o \
  write_shirley_ob.o \
  obf_xml_write.o \


QEMODS = $(BASEMODS)
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a 
MODOBJ = ../../Modules/libqemod.a

TLDEPS= pwlibs

all : tldeps obf_basis.x

obf_basis.x : obf_basis.o libobf.a  $(OBFOBJS) $(PWOBJS) $(QEMODS) $(GWWOBJ)
	$(LD) $(LDFLAGS) -o $@ \
		obf_basis.o libobf.a $(PWOBJS) $(GWWOBJ) $(QEMODS) $(QELIBS) $(LIBMIN)
	- ( cd ../../bin ; ln -fs ../OBF/obf_basis/$@ . )

tldeps :
	if test -n "$(TLDEPS)" ; then   ( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi


libobf.a : $(OBFOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@


clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L libobf.a

include make.depend
# DO NOT DELETE

