include ../make.flags
include modules.flags

all: filter_sequence.x lenz_cat.x ephem_table.x Utest_table.x # convert_order_results.x:

# Executables

LIBDEPS= ../lib/libprop.a ../lib/libmoid.a ../lib/libgauss.a ../lib/libsuit.a

filter_sequence.x: filter_sequence.o option_file.o fildiropn.o spline.o \
	read_interpol_ephem.o filter1.o filter2.o filter3.o \
	filter_output.o $(LIBDEPS)
	$(FC) $(FFLAGS) filter_sequence.o option_file.o fildiropn.o spline.o \
	read_interpol_ephem.o filter1.o filter2.o filter3.o filter_output.o \
	$(LIBDEPS) -o filter_sequence.x

lenz_cat.x: option_file.o lenz_cat.o elem_propag.o $(LIBDEPS)
	$(FC) $(FFLAGS) option_file.o lenz_cat.o elem_propag.o \
	$(LIBDEPS) -o lenz_cat.x

ephem_table.x: ephem_table.o  fildiropn.o $(LIBDEPS)
	$(FC) $(FFLAGS) fildiropn.o  ephem_table.o \
	$(LIBDEPS) -o  ephem_table.x

Utest_table.x: option_file.o fildiropn.o spline.o read_interpol_ephem.o  Utest_table.o $(LIBDEPS)
	$(FC) $(FFLAGS) option_file.o fildiropn.o spline.o read_interpol_ephem.o \
	Utest_table.o $(LIBDEPS) -o Utest_table.x

convert_order_results.x: convert_order_results.o $(LIBDEPS)
	$(FC) $(FFLAGS) convert_order_results.o $(LIBDEPS) -o convert_order_results.x

depend:
	../../bin/mkdep90 -I../include -L../suit -L../propag  *.f *.f90 *.h *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit ../propag  > modules.flags

clean:
	-(rm -f *.x *.o *% *~ *.mod core)

distclean: clean
	../lib/deinst_src


# Dependencies
include make.dep


