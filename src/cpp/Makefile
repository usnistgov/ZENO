# ================================================================
# 
# Disclaimer:  IMPORTANT:  This software was developed at the National
# Institute of Standards and Technology by employees of the Federal
# Government in the course of their official duties.  Pursuant to
# title 17 Section 105 of the United States Code this software is not
# subject to copyright protection and is in the public domain.  This
# is an experimental system.  NIST assumes no responsibility
# whatsoever for its use by other parties, and makes no guarantees,
# expressed or implied, about its quality, reliability, or any other
# characteristic.  We would appreciate acknowledgement if the software
# is used.  This software can be redistributed and/or modified freely
# provided that any derivative works bear some notice that they are
# derived from it, and any modified versions bear some notice that
# they have been modified. 
# 
# ================================================================

# ================================================================
# 
# Author: Derek Juba <derek.juba@nist.gov>
# Date:   Wed Feb 19 11:56:36 2014 EDT
# 
# Time-stamp: <2016-09-28 16:16:36 dcj>
# 
# ================================================================

NANOFLANN_DIR := /home/dcj/Programs/nanoflann-1.1.8
SPRNG_DIR     := /home/dcj/Programs/sprng5

MPI_CXX := mpic++

#MPI_CXXFLAGS := -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi
#MPI_LDFLAGS  := -L/usr/lib/openmpi/lib -lmpi_cxx -lmpi -ldl -lhwloc

GENGETOPT := true #gengetopt

FLEXCXX  := true #flexc++
BISONCXX := true #bisonc++

LEXER   := lexer
GRAMMAR := grammar

PARSERBASE.h := Parserbase.h
PARSE.cc     := parse.cc

SCANNERBASE.h := Scannerbase.h
SCAN.cc       := lex.cc

DEPEND := depend

OPTIONS := options.ggo
CMDLINE := cmdline

CFLAGS   := -Wall -g -O3
CXXFLAGS := $(CFLAGS) -std=c++11

LDFLAGS := -pthread -lz

SCM := SphereCenterModel
RNG := RandomNumber

SOURCES.cc := $(SCAN.cc) $(PARSE.cc) Main.cc Parser.cc Timer.cc ResultsZeno.cc ResultsInterior.cc ResultsCompiler.cc Parameters.cc Units.cc

SOURCES.c := $(CMDLINE).c

############
# RNG Type #
############

#SPRNG RNG
CXXFLAGS   += -I$(SPRNG_DIR)/include
LDFLAGS    += -lsprng -L$(SPRNG_DIR)/lib
SOURCES.cc += $(RNG)/SPRNG.cc

########################
# Sphere Centers Model #
########################

#NanoFLANN Sort SCM
CXXFLAGS   += -I$(NANOFLANN_DIR)/include
SOURCES.cc += $(SCM)/NanoFLANNSort.cc

# ================================================================

OBJECTS.cc := $(SOURCES.cc:.cc=.o)
OBJECTS.c  := $(SOURCES.c:.c=.o)

MPI_OBJECTS.cc := $(SOURCES.cc:.cc=.mpi_o) 

EXECUTABLE := zeno

MPI_EXECUTABLE := zeno-mpi

# ================================================================

.PHONY: all,mpi,clean,distclean,maintainer-clean,check

# ================================================================

all: $(EXECUTABLE)

mpi: $(MPI_EXECUTABLE)

clean:
	$(RM) $(OBJECTS.cc) $(OBJECTS.c) $(MPI_OBJECTS.cc) $(EXECUTABLE) $(MPI_EXECUTABLE)

distclean: clean
	$(RM) $(DEPEND) 

maintainer-clean: distclean
	$(RM) $(CMDLINE).c $(CMDLINE).h $(SCANNERBASE.h) $(SCAN.cc) $(PARSERBASE.h) $(PARSE.cc)

check:
	cd SelfTests && ./run_tests.sh

$(EXECUTABLE): $(OBJECTS.cc) $(OBJECTS.c)
	$(CXX) $(OBJECTS.cc) $(OBJECTS.c) $(LDFLAGS) -o $@

$(MPI_EXECUTABLE): $(MPI_OBJECTS.cc) $(OBJECTS.c)
	$(MPI_CXX) $(MPI_OBJECTS.cc) $(OBJECTS.c) $(LDFLAGS) $(MPI_LDFLAGS) -o $@

$(MPI_OBJECTS.cc): %.mpi_o: %.cc
	$(MPI_CXX) -DUSE_MPI -c $(CXXFLAGS) $(MPI_CXXFLAGS) $< -o $@

$(CMDLINE).c $(CMDLINE).h: $(OPTIONS)
	$(GENGETOPT) --file-name=$(CMDLINE) < $<

$(SCANNERBASE.h) $(SCAN.cc): $(LEXER)
	$(FLEXCXX) $(LEXER)

$(PARSERBASE.h) $(PARSE.cc): $(GRAMMAR)
	$(BISONCXX) $(GRAMMAR)

$(DEPEND): $(SOURCES.cc) $(SOURCES.c)
	$(CXX) -MM -MG $(CXXFLAGS) $^ > $@

# ================================================================

-include $(DEPEND)

# ================================================================

# Local Variables:
# time-stamp-line-limit: 30
# End:
