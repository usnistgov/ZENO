# ================================================================
# 
# This software was developed at the National Institute of
# Standards and Technology by employees of the Federal Government
# in the course of their official duties.  Pursuant to Title 17
# Section 105 of the United States Code, this software is not
# subject to copyright protection and is in the public domain.
# This software is an experimental system.  NIST assumes no
# responsibility whatsoever for its use by other parties, and makes
# no guarantees, expressed or implied, about its quality,
# reliability, or any other characteristic.  We would appreciate
# acknowledgment if the software is used.
#
# ================================================================
# 
# Author: Walid Keyrouz (walid.keyrouz@nist.gov)
# Date:   Tue Jul  2 13:44:54 2013 EDT
# 
# Time-stamp: <2015-03-04 16:48:37 walid>
# 
# ================================================================

#-- From image-stitching/apps/Makefile.in

OS_NAME := $(shell uname -s)
OS_VERS := $(shell uname -r)

ifeq ($(OS_NAME), Darwin)
  F08_C := gfortran-mp-4.8
else
  F08_C := gfortran
endif

# ifeq ($(findstring Darwin, $(OS_NAME)), Darwin)
#   ifeq ($(findstring 12.2.1, $(OS_VERS)), 12.2.1)
#     BOOST_DIR := /usr
#     BOOST_LIB_SUFFIX := -mt
#   else
#     BOOST_DIR := /opt
#     BOOST_LIB_SUFFIX :=
#     GGO_CXXFLAGS += -Wno-write-strings
#   endif
#   BOOST_INCDIR := -I$(BOOST_DIR)/local/include
#   BOOST_LIBDIR := -L$(BOOST_DIR)/local/lib
# else
#   BOOST_INCDIR :=
#   BOOST_LIBDIR :=
# endif

#-- http://fortran90.org/src/faq.html

# Development:

# F08_C := gfortran

F08_PROF_FLAGS := -pg -Wno-omit-frame-pointer

F08_FLAGS := 
F08_FLAGS += -cpp
F08_FLAGS += -Wall -Wextra
F08_FLAGS += -g
F08_FLAGS += -Wsurprising
F08_FLAGS += -Wimplicit-interface
F08_FLAGS += -fPIC -fmax-errors=1
F08_FLAGS += -fcheck=all -fbacktrace
# F08_FLAGS += -Wno-unused-dummy-argument
# F08_FLAGS += -std=f2008 -ffree-form 

F08_FLAGS += $(shell pkg-config --cflags fgsl)

F08_FLAGS += -Wno-compare-reals

# Production:
# -Wall -Wextra -Wimplicit-interface -fPIC -Werror -fmax-errors=1 -O3 -march=native -ffast-math -funroll-loops

%.o: %.f08
	$(F08_C) $(F08_FLAGS) -c $< -o $@

%.pg-o: %.f08
	$(F08_C) $(F08_PG_FLAGS) $(F08_FLAGS) -c $< -o $@

# ----------------------------------------------------------------

LIB_DIR := ../lib
LIBS    += $(shell pkg-config --libs-only-l fgsl)

LDFLAGS += -L$(LIB_DIR)
LDFLAGS += $(shell pkg-config --libs-only-L fgsl)

# ================================================================

EXECS := zeno-v4 zeno-v4-pg

## F08_SRCS_MODS := 
## F08_SRCS_MODS += zeno-main.f08
## F08_SRCS_MODS += math-utils.f08
## F08_SRCS_MODS += stat-utils.f08
## F08_SRCS_MODS += numeric-kinds.f08
## F08_SRCS_MODS += numeric-constants.f08
## F08_SRCS_MODS += zeno-enums.f08
## F08_SRCS_MODS += zeno-parser.f08
## F08_SRCS_MODS += zeno-rng.f08
## F08_SRCS_MODS += zeno-sphere.f08
## F08_SRCS_MODS += zeno-vectors.f08
## F08_SRCS_MODS += zeno-zeerot.f08
## F08_SRCS_MODS += zeno-kernels.f08
## F08_SRCS_MODS += zeno-units.f08
## F08_SRCS_MODS += zeno-tensors.f08
## F08_SRCS_MODS += zeno-toss.f08
## F08_SRCS_MODS += zeno-shape-dist.f08
## F08_SRCS_MODS += zeno-output.f08
## F08_SRCS_MODS += zeno-setup.f08
## F08_SRCS_MODS += zeno-ellipsoid-utils.f08
## F08_SRCS_MODS += zeno-shape-interior-test.f08
## F08_SRCS_MODS += zeno-shape-surf-area.f08
## F08_SRCS_MODS += zeno-shape-volume.f08
## F08_SRCS_MODS += sort-utils.f08
## F08_SRCS_MODS += zeno-shape-surf-points.f08
## F08_SRCS_MODS += zeno-prob-geom.f08
## F08_SRCS_MODS += zeno-integrations.f08
## F08_SRCS_MODS += zeno-debug.f08
## F08_SRCS_MODS += zeno-options.f08
## F08_SRCS_MODS += zeno-warnings.f08
## F08_SRCS_MODS += zeno-codes-data.f08
## F08_SRCS_MODS += zeno-files-data.f08
## F08_SRCS_MODS += getopt.f08
## F08_SRCS_MODS += zeno-eigen.f08
## F08_SRCS_MODS += zeno-version.f08

F08_SRCS_MODS := 

# 0
F08_SRCS_MODS += getopt.f90
F08_SRCS_MODS += numeric-kinds.f90
F08_SRCS_MODS += zeno-codes-data.f90
F08_SRCS_MODS += zeno-enums.f90
F08_SRCS_MODS += zeno-files-data.f90
F08_SRCS_MODS += zeno-version.f90
F08_SRCS_MODS += zeno-warnings.f90

# 1
F08_SRCS_MODS += math-utils.f90
F08_SRCS_MODS += numeric-constants.f90
F08_SRCS_MODS += stat-utils.f90
F08_SRCS_MODS += zeno-debug.f90
F08_SRCS_MODS += zeno-options.f90

# 2
F08_SRCS_MODS += sort-utils.f90
F08_SRCS_MODS += zeno-eigen.f90
F08_SRCS_MODS += zeno-kernels.f90
F08_SRCS_MODS += zeno-rng.f90
F08_SRCS_MODS += zeno-shape-volume.f90
F08_SRCS_MODS += zeno-units.f90
F08_SRCS_MODS += zeno-vectors.f90

# 3
F08_SRCS_MODS += zeno-ellipsoid-utils.f90
F08_SRCS_MODS += zeno-shape-interior-test.f90
F08_SRCS_MODS += zeno-shape-surf-area.f90
F08_SRCS_MODS += zeno-sphere.f90
F08_SRCS_MODS += zeno-tensors.f90
F08_SRCS_MODS += zeno-toss.f90

# 4
F08_SRCS_MODS += zeno-output.f90
F08_SRCS_MODS += zeno-shape-dist.f90
F08_SRCS_MODS += zeno-shape-surf-points.f90
F08_SRCS_MODS += zeno-zeerot.f90

# 5
F08_SRCS_MODS += zeno-integrations.f90
F08_SRCS_MODS += zeno-parser.f90
F08_SRCS_MODS += zeno-prob-geom.f90
F08_SRCS_MODS += zeno-setup.f90

# 6
F08_SRCS_MODS += zeno-main.f90

F08_OBJS := $(F08_SRCS_MODS:.f90=.o)

FILES_GEN :=
FILES_GEN += .deps_list.mk .mods_list.mk

SYS_MODS :=
SYS_MODS += fgsl

# ================================================================

.phony: all clean reallyclean modules mods.list deps.list tags

# ----------------------------------------------------------------

all: $(EXECS)

# all: zeno-v4

# ----------------------------------------------------------------

-include .mods_list.mk
-include .deps_list.mk

# ----------------------------------------------------------------

modules: .mods_list.mk .deps_list.mk $(F08_SRCS_MODS:.f90=.o)

# ----------------------------------------------------------------

zeno-v4: $(F08_SRCS_MODS:.f08=.o)
	$(F08_C) -o $@ $^ $(LDFLAGS) $(LIBS)

zeno-v4-pg: $(F08_SRCS_MODS) zeno-v4
	$(F08_C) $(F08_PROF_FLAGS) -o $@ $(F08_FLAGS) $(F08_SRCS_MODS) $(LDFLAGS) $(LIBS)

# ----------------------------------------------------------------

mods.list : .mods_list.mk

.mods_list.mk : ${F08_SRCS_MODS}
	gawk -f ./mkdepend.awk -v EMIT_MOD_LIST=1 $^ > .mods_list.mk
#	./f08_deps.py -m modules-list -x fgsl $^ > $@

deps.list : .deps_list.mk

.deps_list.mk : ${F08_SRCS_MODS}
	gawk -f ./mkdepend.awk -v EXCLUDE_MODS="$(SYS_MODS)" $^ > .deps_list.mk
#	./f08_deps.py -x fgsl $^ > $@

# ----------------------------------------------------------------

tags: TAGS

TAGS: $(F08_SRCS_MODS)
	etags -l fortran $^

# ----------------------------------------------------------------

clean:
	- $(RM) $(EXECS) $(wildcard *.o) $(wildcard *.mod)

reallyclean: clean
	- $(RM) $(FILES_GEN) TAGS

# ================================================================

# Local Variables:
# time-stamp-line-limit: 30
# mode: GNUmakefile
# End:
