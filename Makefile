#Directories
PWD    = $(shell pwd)
OBJDIR = ./obj
SRCDIR = ./src/main
#LIBDIR = $(PWD)/UTILS/lib
#IDIR   = $(PWD)/UTILS/include
#EVALLIB=/home/lei/ASDF/evalresp/.libs

CFLAGS = -check all,noarg_temp_created

# new version with default sac libraries
TAULIBDIR=$(PWD)/ttimes_mod
SACLIBDIR = ${SACHOME}/lib
LIBS = -lsacio -lsac -lm -ltau -lasdf
 #LIB = -L/opt/seismo/lib -lDRWFiles -lf90recipes -lDSacio -lDSacLib -lSacTools -lm

#all_obj = $(shell find . -name obj/*.*)

ADIOS_FLIB=$(shell adios_config -lf)
ADIOS_INC=$(shell adios_config -cf)

ASDFLIBDIR=$(ASDFHOME)/lib
ASDFINCDIR=$(ASDFHOME)/include

## set ADIOS_DIR here or before doing make
#override ADIOS_DIR:=/home/lei/bin/adios-1.5.0
#override ADIOS_INC:=`${ADIOS_DIR}/bin/adios_config -c -f`
#override ADIOS_FLIB:=`${ADIOS_DIR}/bin/adios_config -l -f`

############################
#compiler option
#OPT = -I${SHARED}
#OPT = -std03
FC = ifort
CC = icc
MPIFC = mpif90
MPICC = mpicc
CFLAGS= -g -O0 -check noarg_temp_created

_OBJ = main_subs.o main.o 

OBJ = $(patsubst %, ${OBJDIR}/%, $(_OBJ))

##########################################################
PROG = Measure_adj 
default: MK_OBJDIR ${PROG}

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	  $(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) -I$(ASDFINCDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.f
	  $(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) -I$(ASDFINCDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	  $(MPICC) -c -o $@ $< 

#include shared/Makefile FLEXWIN/Makefile measure_adj/Makefile

MK_OBJDIR:
	mkdir -p $(OBJDIR)

make_shared:
	cd src/shared; make

make_ma:
	cd src/measure_adj; make

all_obj = $(wildcard $(OBJDIR)/*.o)

${PROG}: make_shared make_ma $(OBJ)
	${MPIFC} ${CFLAGS} -o $@ $(all_obj) \
		-L${TAULIBDIR} -L${SACLIBDIR} -L$(ASDFLIBDIR) ${LIBS} ${ADIOS_FLIB}


.PHONY:clean print_var cleanall

print_var:
	@echo $(OBJ)
	@echo $(SRCDIR)
	@echo $(all_obj)

clean:
	rm -f  ${LIB_ALL} ${PROG} *.o *.mod *.a $(OBJDIR)/*

cleanall:
	rm -f  iasp91.*
	cd ${TAULIBDUR} ; make -f make_gfortran clean

