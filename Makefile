#Directories
PWD    = $(shell pwd)
OBJDIR = ./obj
SRCDIR = ./src/main
#LIBDIR = $(PWD)/UTILS/lib
#IDIR   = $(PWD)/UTILS/include
#EVALLIB=/home/lei/ASDF/evalresp/.libs

CFLAGS = -check

# new version with default sac libraries
TAULIBDIR=$(PWD)/ttimes_mod
SACLIBDIR = ${SACHOME}/lib
LIBS = -lsacio -lsac -lm -ltau
 #LIB = -L/opt/seismo/lib -lDRWFiles -lf90recipes -lDSacio -lDSacLib -lSacTools -lm

#all_obj = $(shell find . -name obj/*.*)

## set ADIOS_DIR here or before doing make
override ADIOS_DIR:=/home/lei/bin/adios-1.5.0
override ADIOS_INC:=`${ADIOS_DIR}/bin/adios_config -c -f`
override ADIOS_FLIB:=`${ADIOS_DIR}/bin/adios_config -l -f`

############################
#compiler option
#OPT = -I${SHARED}
#OPT = -std03
FC = ifort
CC = icc
MPIFC = mpif90
MPICC = mpicc
CFLAGS= -g -O0 -check all

_OBJ = mpi_weighting_subs.o rotate_subs.o var_main_mod.o main_subs.o main.o 

OBJ = $(patsubst %, ${OBJDIR}/%, $(_OBJ))

## set ADIOS_DIR here or before doing make
override ADIOS_DIR:=/home/lei/bin/adios-1.5.0
override ADIOS_INC:=` ${ADIOS_DIR}/bin/adios_config -c -f`
override ADIOS_FLIB:=`${ADIOS_DIR}/bin/adios_config -l -f`

##########################################################
PROG = TEST
default: ${PROG}

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	  $(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) $(ADIOS_INC)

$(OBJDIR)/%.o: $(SRCDIR)/%.f
	  $(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) $(ADIOS_INC)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	  $(MPICC) -c -o $@ $< 

#include shared/Makefile FLEXWIN/Makefile measure_adj/Makefile

make_asdf:
	cd src/asdf_util; make

make_shared:
	cd src/shared; make

make_ma:
	cd src/measure_adj; make

all_obj = $(wildcard $(OBJDIR)/*.o)

${PROG}: make_asdf make_shared make_ma $(OBJ)
	${MPIFC} ${CFLAGS} -o $@ $(all_obj) -I$(OBJDIR) \
		-L${TAULIBDIR} -L${SACLIBDIR} ${LIBS} ${ADIOS_FLIB}


.PHONY:clean print_var cleanall

print_var:
	@echo $(OBJ)
	@echo $(SRCDIR)
	@echo $(all_obj)

clean:
	rm -f  ${LIB_ALL} ${PROGS} *.o *.mod *.a $(OBJDIR)/*

cleanall:
	rm -f  iasp91.*
	cd ${TAULIBDUR} ; make -f make_gfortran clean

