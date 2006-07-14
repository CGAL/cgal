# Copyright (c) 1997-2006  Utrecht University (The Netherlands),
# ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
# INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
# (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
# and Tel-Aviv University (Israel).  All rights reserved.
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
# Author(s)     : Radu Ursu
#                 Sylvain Pion

# This is the makefile for compiling static and shared libraries in CGAL.

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
assertions$(OBJ_EXT) \
Bbox_2_intersections$(OBJ_EXT) \
Bbox_3_intersections$(OBJ_EXT) \
cgal_logo$(OBJ_EXT) \
Color$(OBJ_EXT) \
File_header_extended_OFF$(OBJ_EXT) \
File_header_OFF$(OBJ_EXT) \
File_scanner_OFF$(OBJ_EXT) \
File_writer_inventor$(OBJ_EXT) \
File_writer_OFF$(OBJ_EXT) \
File_writer_VRML_2$(OBJ_EXT) \
File_writer_wavefront$(OBJ_EXT) \
Geomview_stream$(OBJ_EXT) \
Interval_arithmetic$(OBJ_EXT) \
io$(OBJ_EXT) \
JAMA_numeric_solver$(OBJ_EXT) \
KDS_Log$(OBJ_EXT) \
kernel$(OBJ_EXT) \
MP_Float$(OBJ_EXT) \
numeric_solvers_support$(OBJ_EXT) \
NefPolynomial$(OBJ_EXT) \
Random$(OBJ_EXT) \
Real_timer$(OBJ_EXT) \
Timer$(OBJ_EXT) \
Turkowski_numeric_solver$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

shared_lib: $(CGAL_OBJECTS)
	$(CGAL_SHARED_LIB_CREATE)$(CGAL_SHARED_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_SHARED_LIB_LDFLAGS)
	mv $(CGAL_SHARED_LIB) '$(CGAL_LIB_DESTINATION)'
	rm $(CGAL_OBJECTS)

static_lib: static_lib_no_install
	mv $(CGAL_LIB) $(CGAL_LIB_DESTINATION)

static_lib_no_install: $(CGAL_OBJECTS)
	$(CGAL_LIB_CREATE)$(CGAL_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_LIB_LDFLAGS)
	$(RANLIB) $(CGAL_LIB)
	rm $(CGAL_OBJECTS)

clean:
	rm -f $(CGAL_LIB) $(CGAL_SHARED_LIB) $(CGAL_OBJECTS)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.SUFFIXES: .cpp $(SUFFIXES)

.cpp$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS)  -c $<

