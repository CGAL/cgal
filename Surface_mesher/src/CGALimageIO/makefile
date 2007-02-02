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
  ImageIO$(OBJ_EXT) \
  analyze$(OBJ_EXT) \
  bmp$(OBJ_EXT) \
  bmpendian$(OBJ_EXT) \
  bmpread$(OBJ_EXT) \
  convert$(OBJ_EXT) \
  gif$(OBJ_EXT) \
  gis$(OBJ_EXT) \
  inr$(OBJ_EXT) \
  iris$(OBJ_EXT) \
  mincio$(OBJ_EXT) \
  pnm$(OBJ_EXT) \
  recbuffer$(OBJ_EXT) \
  recline$(OBJ_EXT) \
  reech4x4$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

shared_lib: $(CGAL_OBJECTS)
	$(CGAL_IMAGEIO_SHARED_LIB_CREATE)$(CGAL_IMAGEIO_SHARED_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_SHARED_LIB_LDFLAGS)
	mv $(CGAL_IMAGEIO_SHARED_LIB) '$(CGAL_LIB_DESTINATION)'
	rm $(CGAL_OBJECTS)

static_lib: static_lib_no_install
	mv $(CGAL_IMAGEIO_LIB) $(CGAL_LIB_DESTINATION)

static_lib_no_install: $(CGAL_OBJECTS)
	$(CGAL_LIB_CREATE)$(CGAL_IMAGEIO_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_LIB_LDFLAGS)
	$(RANLIB) $(CGAL_IMAGEIO_LIB)
	rm $(CGAL_OBJECTS)

clean:
	rm -f $(CGAL_IMAGEIO_LIB) $(CGAL_IMAGEIO_SHARED_LIB) $(CGAL_OBJECTS)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.SUFFIXES: .cpp $(SUFFIXES)

.cpp$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS)  -c $<

.c$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS)  -c $<

