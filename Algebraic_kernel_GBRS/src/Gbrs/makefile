# Copyright (c) 2006 Inria Lorraine (France). All rights reserved.
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

# This is the makefile for compiling the CGAL GB/RS support.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# The following line should better be passed as command line argument.
# CGAL_MAKEFILE = SHOULD_BE_SET_BY_INSTALL_CGAL
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

# CXXFLAGS needs to be passed as command line argument
# CXXFLAGS = $(CGAL_SHARED_LIB_CXXFLAGS)

#---------------------------------------------------------------------#
#                    object files
#---------------------------------------------------------------------#

CGAL_OBJECTS = \
Gbrs_algebraic_1$(OBJ_EXT) \
Gbrs_polynomial_1$(OBJ_EXT) \
Gbrs_solve_1$(OBJ_EXT) \
Gbrs_polynomial_2$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: outside

outside: CGAL_GBRS_LIB=libCGALGbrs.a
outside: CXXFLAGS=$(CGAL_CXXFLAGS)
outside: OUTSIDEPARS=-I$(PWD)/../../include -I$(MPFI_INCL) -I$(RS_INCL)
outside: static_lib_no_install
	mv $(CGAL_GBRS_LIB) $(CGAL_LIB_DESTINATION)

static_lib: static_lib_no_install
	mv $(CGAL_GBRS_LIB) $(CGAL_LIB_DESTINATION)

static_lib_no_install: $(CGAL_OBJECTS)
	$(CGAL_LIB_CREATE)$(CGAL_GBRS_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_LIB_LDFLAGS)
	$(RANLIB) $(CGAL_GBRS_LIB)
	rm $(CGAL_OBJECTS)

clean:
	rm -f $(CGAL_GBRS_LIB) $(CGAL_OBJECTS)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.SUFFIXES: .cpp $(SUFFIXES)

.cpp$(OBJ_EXT):
	$(CGAL_CXX)  $(OUTSIDEPARS)  $(CXXFLAGS)  -c $<

