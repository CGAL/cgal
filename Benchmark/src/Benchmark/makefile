# Copyright (c) 2006 Tel-Aviv University (Israel).
# All rights reserved.
#
# This file is part of CGAL (www.cgal.org); you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; version 2.1 of the License.
# See the file LICENSE.LGPL distributed with CGAL.
#
# Licensees holding a valid commercial license may use this file in
# accordance with the commercial license agreement provided with the software.
#
# This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
# WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# $URL$
# $Id$
#
# Author(s): Efi Fogel         <efif@post.tau.ac.il>

# This is the makefile for compiling the Benchmark support library libCGALBenchmark.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = SHOULD_BE_SET_BY_COMMAND_LINE

include $(CGAL_MAKEFILE)

# defaut rule
all: lib

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

# CXXFLAGS must be provided on the make command-line.
# CXXFLAGS = -I"../../include" $(CGAL_WINLIB_CXXFLAGS) -DCGAL_USE_BENCHMARK

#---------------------------------------------------------------------#
# Object files
#---------------------------------------------------------------------#

CGAL_OBJECTS = \
Benchmark$(OBJ_EXT) \
benchmark_parser$(OBJ_EXT) \
benchmark_lexer$(OBJ_EXT) \
Benchmark_visitor$(OBJ_EXT) \
Option_parser$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

# TODO : shared_lib !

static_lib: lib

lib: lib_no_install
	mv $(CGAL_BENCHMARK_LIB) $(CGAL_LIB_DESTINATION)

lib_no_install: $(CGAL_OBJECTS) 
	$(CGAL_LIB_CREATE)$(CGAL_BENCHMARK_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_LIB_LDFLAGS)
	$(RANLIB) $(CGAL_BENCHMARK_LIB)
	rm $(CGAL_OBJECTS) 

clean:
	rm -f $(CGAL_BENCHMARK_LIB) $(CGAL_OBJECTS)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

%$(OBJ_EXT): %.cpp
	$(CGAL_CXX) $(CXXFLAGS) -c $<

.SUFFIXES: .cpp $(SUFFIXES)

.cpp$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS) -c $<

FLEX ?=flex
BISON ?=bison

# Apparently, old versions of bison (e.g., 1.28) name the generated definion
# header file <base>.cpp.h. The file must be renamed to <base>.hpp 
BISON_VERSION_CMD =expr match "$$(bison --version)" '.*\([1-9]\.[0-9]*\)'
BISON_VERSION =$(shell $(BISON_VERSION_CMD))
OLD_BISON_VERSION_CMD =expr "$(BISON_VERSION)" \<= 1.28
OLD_BISON_VERSION =$(shell $(OLD_BISON_VERSION_CMD))

%.cpp %.hpp : %.y
	$(BISON) -d $< -o $*.cpp
ifeq ($(OLD_BISON_VERSION), 1)
	mv $*.cpp.h $*.hpp
endif

%.cpp : %.l
	$(FLEX) -8 -o$@ $<
