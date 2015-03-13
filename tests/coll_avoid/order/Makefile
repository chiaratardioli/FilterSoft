COMPILER = ifort
OPT = -g -CB -CA -CU -traceback -debug all -O3 # -warn all -e90 -ftrapuv 
EXECUTABLE = order

all: clearScreen compile exe clean plot convert

bau: clearScreen exe clean plot convert

clearScreen:
	clear

compile: modules
	$(COMPILER) $(OPT) *.o -o $(EXECUTABLE)

modules:
	$(COMPILER) $(OPT) -c parameters.f90 geom.f90 pol_zeros.f90 conv_mslib.f90 jplephem.f geopotential.f90 julian.f90 polynom.f90 tbp.f90 radPress.f90 system.f90 integrator.f90 main.f90

exe: 
	./$(EXECUTABLE)

plot: 
	gnuplot plots.gnu -

convert:
	ps2pdf -dEPSCrop cart.eps
	ps2pdf -dEPSCrop kepl.eps
	ps2pdf -dEPSCrop cart3D.eps	

clean:
	rm -rf *.o *.mod *~ *.dat \#*\#






