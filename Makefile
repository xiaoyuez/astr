# This makefile is used to compile ASTR code.
# The compiler: gfortran compiler
#
#FCFLAGS=  -Wuse-without-only -g
#FC=mpif90
FC=h5pfc

SRCDIR = src
OBJDIR = obj
BINDIR = bin
# CTRDIR =  /home/fangjian/opt/cantera_2.5.1

FCFLAGS= -O3 -fbounds-check

# OPTIONS1 = -fcheck=all
OPTIONS2 = -J $(OBJDIR)
OPTIONS3 = -DHDF5
# OPTIONS4 = -DCOMB -I$(CTRDIR)/include/cantera 
# OMP = -fopenmp

EXE=astr

# LIBS= -lz -lm -L$(CTRDIR)/lib -lcantera_fortran -lcantera -lstdc++ -pthread
LIBS= -lz -lm 

TARGET = $(BINDIR)/$(EXE)

VPATH = $(SRCDIR):$(OBJDIR)

srs=  fdnn.F90 singleton.F90 commtype.F90 stlaio.F90 constdef.F90 tecio.F90 vtkio.F90     \
      interp.F90 commvar.F90 thermchem.F90 commarray.F90 fludyna.F90 parallel.F90         \
      hdf5io.F90 cmdefne.F90 commfunc.F90 commcal.F90 models.F90 statistic.F90 bc.F90     \
      readwrite.F90 geom.F90 ibmethod.F90 gridgeneration.F90 riemann.F90 solver.F90       \
      pp.F90 initialisation.F90 mainloop.F90 astr.F90
OBJS=$(srs:.F90=.o)

%.o:%.F90
	@mkdir -p $(OBJDIR) 
	$(FC) $(FCFLAGS) $(INCL) $(OPTIONS1) $(OPTIONS2) $(OPTIONS3) $(OPTIONS4) $(OMP) -c -o $(OBJDIR)/$@  $<

default: $(OBJS)
	@mkdir -p $(BINDIR)
	$(FC) $(FCFLAGS) -o $(TARGET) $(OBJDIR)/*.o $(LIBS) $(INCL) $(OMP)

clean:
	rm -fv $(OBJDIR)/*.o $(OBJDIR)/*.mod $(TARGET) $(SRCDIR)/*.mod