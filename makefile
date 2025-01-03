CPL?=g++
OPTLVL=-O3
DEBUGCPP=-g
CPLOPT=-std=c++11 -Wall $(DEBUGCPP) $(OPTLVL) -fPIC
FCC?=gfortran
PYTHON?=python3

LINKER?=g++

ifeq ($(FCC), ifort)
    lflag_prepro ?= -fpp
    lflags_precise = -fp-model precise -fimf-arch-consistency=true $(lflag_prepro) -fPIC
    lflags_debug = -debug -save-temps -fpic -heap-arrays -O0 -g -traceback -check all -fpe0 -fp-stack-check $(lflags_precise)
    lflags_fast = -O2 $(lflag_prepro)
    lflags_profiling = $(lflags_fast) -profile-functions -profile-loops=all -profile-loops-report=2 $(lflag_prepro)
    LINKOPT?=-lifort
else
ifeq ($(FCC), gfortran)
    lflag_prepro ?= -cpp
    lflags_precise = -g -Wall -fcheck=all -ffpe-trap=zero,invalid,overflow $(lflag_prepro) -fPIC
    lflags_debug = $(lflags_precise) -fbacktrace -ffree-line-length-0
    lflags_fast = -O3 $(lflag_prepro)
    lflags_profiling = $(lflags_fast) -pg $(lflag_prepro)
    LINKOPT?=-L"$(shell dirname `gfortran --print-file-name libgfortran.a`)" -lgfortran
else
    $(error UnknownCompiler: $(FCC))
endif
endif

FCCOPT?=-c $(lflags_debug)
FCCOPTLEGACY?=$(FCCOPT) -std=legacy


exe?=re

all: $(exe)

OBJS=main.o rate_equation_lsode.o ode_wrapper.o types.o constants.o opkdmain.o opkda1.o opkda2.o calculate_reaction_rate.o utils.o

OBJSWRAPPER=rate_equation_lsode.o ode_wrapper.o types.o constants.o calculate_reaction_rate.o utils.o opkdmain.o opkda1.o opkda2.o
CONSTSOBJSWRAPPER=constants.o

LIB_DIR=.
CHEMPLLIB=$(LIB_DIR)/libchempl.a
MYCONSTSLIB=$(LIB_DIR)/libmyconsts.a

chempl: setup.py chempl.pyx myconsts.pyx myconsts.pxd $(CHEMPLLIB) $(MYCONSTSLIB)
	$(PYTHON) setup.py build_ext --inplace

$(CHEMPLLIB): $(OBJSWRAPPER)
	ar rcs $(CHEMPLLIB) $(OBJSWRAPPER)

$(MYCONSTSLIB): $(CONSTSOBJSWRAPPER)
	ar rcs $(MYCONSTSLIB) $(CONSTSOBJSWRAPPER)

$(exe): makefile $(OBJS)
	$(LINKER) $(OBJS) $(LINKOPT) -o $(exe)

types.o: types.hpp types.cpp
	$(CPL) $(CPLOPT) -c types.cpp

constants.o: constants.cpp
	$(CPL) $(CPLOPT) -c constants.cpp

main.o: main.cpp rate_equation_lsode.hpp types.hpp calculate_reaction_rate.hpp
	$(CPL) $(CPLOPT) -c main.cpp

calculate_reaction_rate.o: calculate_reaction_rate.cpp calculate_reaction_rate.hpp
	$(CPL) $(CPLOPT) -c calculate_reaction_rate.cpp

rate_equation_lsode.o: rate_equation_lsode.cpp calculate_reaction_rate.hpp
	$(CPL) $(CPLOPT) -c rate_equation_lsode.cpp

utils.o: utils.cpp
	$(CPL) $(CPLOPT) -c utils.cpp

ode_wrapper.o: ode_wrapper.f90
	$(FCC) $(FCCOPT) ode_wrapper.f90

opkdmain.o: opkdmain.f opkda1.f opkda2.f
	$(FCC) $(FCCOPTLEGACY) opkd*.f

clean:
	rm *.o *.gch *.a *.so
