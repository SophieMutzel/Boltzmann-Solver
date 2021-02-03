# Makefile for Boltzmann code, adapted from pseudospectators/FLUSI and pseudospectators/UP2D
# Non-module Fortran files to be compiled:
FFILES =
# Object and module directory:
OBJDIR = OBJ
OBJS := $(FFILES:%.f90=$(OBJDIR)/%.o)

# Files that create modules:
MFILES = module_precision.f90 module_params.f90 module_read_write.f90 module_utils.f90 \
				module_cosmo.f90 module_xsecs.f90 module_rhs.f90 dvode.f90
MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)

# Source code directories (colon-separated):
VPATH = src

# Set the default compiler if it's not already set
ifndef $(FC)
FC = mpif90
endif

#-------------------------------------------------------------------------------
# COMPILER-DEPENDEND PART
#-------------------------------------------------------------------------------
# GNU compiler
#-------------------------------------------------------------------------------
ifeq ($(shell $(FC) --version 2>&1 | head -n 1 | head -c 3),GNU)
# Specify directory for compiled modules:
FFLAGS += -J$(OBJDIR) # specify directory for modules.
FFLAGS += -O3 -ffree-line-length-none
PPFLAG= -cpp #preprocessor flag
#LDFLAGS = -llapack
# Debug flags for gfortran:
FFLAGS += -Wuninitialized -fimplicit-none -fbounds-check -g -ggdb -pedantic
FFLAGS += -Wall -Wextra -Wconversion -g3 -fbacktrace -ffpe-trap=zero,invalid -finit-real=nan -finit-integer=-99999
FFLAGS += -Wno-unused-variable -Wno-unused-parameter -Wno-unused-dummy-argument -Wno-unused-function
endif

#-------------------------------------------------------------------------------
# Intel compiler
#-------------------------------------------------------------------------------
mpif90:=$(shell $(FC) --version | head -c 5)
ifeq ($(mpif90),ifort)
PPFLAG= -fpp
FFLAGS = -FR -O3 -warn all,nounused -traceback -check bounds -debug all -check all,noarg_temp_created
FFLAGS += -module $(OBJDIR) # specify directory for modules.
#LDFLAGS = -L/usr/X11/lib/ -lX11 #-L/usr/lib64/lapack -llapack
endif

all: directories boltzmann

# Compile main programs, with dependencies.
boltzmann: main.f90 $(MOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ #$(LDFLAGS)

# Compile modules (module dependency must be specified by hand in
# Fortran). Objects are specified in MOBJS (module objects).
$(OBJDIR)/module_precision.o: module_precision.f90
	$(FC) $(FFLAGS) -c -o $@ $< #$(LDFLAGS)

$(OBJDIR)/module_read_write.o: module_read_write.f90 $(OBJDIR)/module_precision.o
	$(FC) $(FFLAGS) -c -o $@ $< #$(LDFLAGS)

$(OBJDIR)/module_utils.o: module_utils.f90 $(OBJDIR)/module_precision.o \
	quadpack.f90 interpolation.f90
	$(FC) $(FFLAGS) -c -o $@ $< #$(LDFLAGS)

$(OBJDIR)/module_params.o: module_params.f90 $(OBJDIR)/module_precision.o $(OBJDIR)/module_read_write.o $(OBJDIR)/module_utils.o \
	ini_cons_to_params.f90 allocate_couplings.f90
	$(FC) $(FFLAGS) -c -o $@ $< #$(LDFLAGS)

$(OBJDIR)/module_cosmo.o: module_cosmo.f90 $(OBJDIR)/module_precision.o $(OBJDIR)/module_params.o $(OBJDIR)/module_utils.o \
	initial_conditions.f90
	$(FC) $(FFLAGS) -c -o $@ $< #$(LDFLAGS)

$(OBJDIR)/module_xsecs.o: module_xsecs.f90 $(OBJDIR)/module_precision.o $(OBJDIR)/module_utils.o
	$(FC) $(FFLAGS) -c -o $@ $< #$(LDFLAGS)

$(OBJDIR)/module_rhs.o: module_rhs.f90 $(OBJDIR)/module_precision.o $(OBJDIR)/module_params.o \
	$(OBJDIR)/module_utils.o $(OBJDIR)/module_xsecs.o $(OBJDIR)/module_cosmo.o \
	rhs_boltzmann.f90 RK4.f90 rhs_region3a2.f90 region3a_log.f90 rhs_contributions.f90 region3aeq.f90
	$(FC) $(FFLAGS) -c -o $@ $< #$(LDFLAGS)

$(OBJDIR)/dvode.o: dvode.f90  $(OBJDIR)/module_params.o $(OBJDIR)/module_utils.o

# Compile remaining objects from Fortran files.
$(OBJDIR)/%.o: %.f90 $(MOBJS)
	$(FC) $(FFLAGS) -c -o $@ $< #$(LDFLAGS)

clean:
	rm -rf $(PROGRAMS) $(OBJDIR) a.out boltzmann

# If the object directory doesn't exist, create it.
.PHONY: directories

directories: ${OBJDIR}

${OBJDIR}:
	mkdir -p ${OBJDIR}
