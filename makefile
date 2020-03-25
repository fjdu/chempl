CPL?=g++
OPTLVL=-O3
DEBUGCPP=-g
CPLOPT=-std=c++11 -Wall $(DEBUGCPP) $(OPTLVL)
FCC?=gfortran

#LINKER?=gfortran
#LINKOPT=-lc++ -Wall $(OPTLVL)
LINKER?=g++
LINKOPT?=-lgfortran -L/usr/local/Cellar/gcc/9.2.0_1/lib/gcc/9/

ifeq ($(FCC), ifort)
    lflag_prepro ?= -fpp
    lflags_precise = -fp-model precise -fimf-arch-consistency=true $(lflag_prepro)
    lflags_debug = -debug -save-temps -fpic -heap-arrays -O0 -g -traceback -check all -fpe0 -fp-stack-check $(lflags_precise)
    lflags_fast = -O2 $(lflag_prepro)
    lflags_profiling = $(lflags_fast) -profile-functions -profile-loops=all -profile-loops-report=2 $(lflag_prepro)
else
ifeq ($(FCC), gfortran)
    lflag_prepro ?= -cpp
    lflags_precise = -g -Wall -fcheck=all -ffpe-trap=zero,invalid,overflow $(lflag_prepro)
    lflags_debug = $(lflags_precise) -fbacktrace -ffree-line-length-0
    lflags_fast = -O3 $(lflag_prepro)
    lflags_profiling = $(lflags_fast) -pg $(lflag_prepro)
else
    $(error UnknownCompiler: $(FCC))
endif
endif

FCCOPT?=-c $(lflags_debug)


exe?=re

all: $(exe)

OBJS=main.o rate_equation_lsode.o ode_wrapper.o types.o constants.o opkdmain.o opkda1.o opkda2.o logistics.o calculate_reaction_rate.o utils.o

OBJSWRAPPER=rate_equation_lsode.o ode_wrapper.o types.o constants.o logistics.o calculate_reaction_rate.o utils.o opkdmain.o opkda1.o opkda2.o

LIB_DIR=.
CHEMPLAYLIB=$(LIB_DIR)/libchemplay.a

chemplay: setup.py chemplay.pyx myconsts.pyx myconsts.pxd $(CHEMPLAYLIB)
	python3 setup.py build_ext --inplace

$(CHEMPLAYLIB): $(OBJSWRAPPER)
	ar rcs $(CHEMPLAYLIB) $(OBJSWRAPPER)

$(exe): $(OBJS)
	$(LINKER) $(LINKOPT) -o $(exe) $(OBJS)

types.o: types.hpp types.cpp
	$(CPL) $(CPLOPT) -c types.cpp

constants.o: constants.cpp
	$(CPL) $(CPLOPT) -c constants.cpp

main.o: main.cpp rate_equation_lsode.hpp logistics.hpp types.hpp calculate_reaction_rate.hpp
	$(CPL) $(CPLOPT) -c main.cpp

logistics.o: logistics.cpp
	$(CPL) $(CPLOPT) -c logistics.cpp

calculate_reaction_rate.o: calculate_reaction_rate.cpp calculate_reaction_rate.hpp
	$(CPL) $(CPLOPT) -c calculate_reaction_rate.cpp

rate_equation_lsode.o: rate_equation_lsode.cpp calculate_reaction_rate.hpp
	$(CPL) $(CPLOPT) -c rate_equation_lsode.cpp

utils.o: utils.cpp
	$(CPL) $(CPLOPT) -c utils.cpp

ode_wrapper.o: ode_wrapper.f90
	$(FCC) $(FCCOPT) ode_wrapper.f90

opkdmain.o: opkdmain.f opkda1.f opkda2.f
	$(FCC) $(FCCOPT) opkd*.f

clean:
	rm *.o *.gch *.a *.so
