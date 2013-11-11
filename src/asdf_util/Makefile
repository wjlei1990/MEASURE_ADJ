#Directories
PWD    = $(shell pwd)
OBJDIR = ../../obj
SRCDIR = $(PWD)
#LIBDIR = $(PWD)/UTILS/lib
#IDIR   = $(PWD)/UTILS/include
#EVALLIB=/home/lei/ASDF/evalresp/.libs

CFLAGS = -check

# new version with default sac libraries
SACLIBDIR = ${SACHOME}/lib
LIBS = -lsacio -lsac -lm
 #LIB = -L/opt/seismo/lib -lDRWFiles -lf90recipes -lDSacio -lDSacLib -lSacTools -lm


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
CFLAGS= -g -check all

_OBJ = adios_helpers_definitions.o adios_helpers_writers.o \
			 adios_helpers.o asdf_data.o asdf_subs.o \

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
.F.o:
	$(OBJDIR)/%.o: $(SRCDIR)/%.F
	  $(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) $(ADIOS_INC)

#$(OBJDIR)/%.o: $(SRCDIR)/%.c
#	  $(MPICC) -c -o $@ $< 

${PROG}: $(OBJ)


.PHONY:clean print_var cleanall

print_var:
	@echo $(OBJ)
	@echo $(SRCDIR)

clean:
	rm -f  ${LIB_ALL} ${PROGS} *.o *.mod *.a $(OBJDIR)/*

cleanall:
	rm -f  iasp91.*
	cd ${TAULIBDUR} ; make -f make_gfortran clean

