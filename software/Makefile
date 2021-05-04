########################################################################
#                                                                      #
# This Makefile generates the executables for NORMA and URO            #
#                                                                      #
########################################################################

PROF      = 

# Options for Linux ####################################################

# GNU g77 compiler
# -static option to generate code independant of dynamically linked libraries
# -m32 option to generate portable 32 bit code 
F77       = g77
F77FLAGS  = -funroll-all-loops -fno-f2c -O3 -static -m32

# Intel FORTRAN compiler : works - but ...
# .. option -f77rtl ( force FORTRAN77 behaviour ) required for NAMELIST reading
# .. when using an Intel-based system, you may add option -fast
# .. to compile NORMA.f using ifort you must delete the declaration of
#    "rand" since the ifort compiler treats RAND as an intrinsic function
# .. there will be no large speed-up of NORMA since most time is spent 
#    in the URO subroutines which are precompiled in 32bit mode under g77
# F77       = ifort
# F77FLAGS  = -f77rtl -O3

# The gfortran compiler does not work with diagrtb (I/O format problems)
# F77       = gfortran
# F77FLAGS  = 

LOADFLGS  = $(F77FLAGS)
ARCH      = LINUX


# Generate executables #################################################
ALL	= diagstd diagrtb proj_modes_bin pdbmat NORMA.exe

all: $(ALL)

proj_modes_bin: proj_modes_bin.f
	$(F77) $(LOADFLGS) -o $@ proj_modes_bin.f

diagstd: diagstd.f
	$(F77) $(LOADFLGS) -o $@ diagstd.f

diagrtb: diagrtb.f
	$(F77) $(LOADFLGS) -o $@ diagrtb.f

pdbmat: pdbmat.f
	$(F77) $(LOADFLGS) -o $@ pdbmat.f

NORMA.exe: NORMA.f
	$(F77) $(LOADFLGS) -o $@ NORMA.f

########################################################################


.f.o:
	$(F77) $(F77FLAGS) -c $?

clean:
	rm -f *.o core func.in func.log func.out *.log 1hhp.eigenfacs 1hhp.proj_modes_bin.IO 1hhp.NMApert.pdb $(ALL)

