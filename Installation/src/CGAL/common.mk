# Copyright (c) 2007  INRIA Sophia-Antipolis (France).
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
# Author(s)     : Laurent Rineau

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
#                    common.mk variables
#---------------------------------------------------------------------#

COMMON_SHARED_LIB = $(SHARED_LIB)$(SOVERSION)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

shared_lib: shared_lib_no_install
	mv $(COMMON_SHARED_LIB) '$(CGAL_LIB_DESTINATION)'


shared_lib_no_install: $(OBJECTS)
	$(CGAL_SHARED_LIB_CREATE)$(COMMON_SHARED_LIB) $(CGAL_SHARED_LIB_SONAME) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_SHARED_LIB_LDFLAGS) $(SHARED_LIB_ADDITIONNAL_LDFLAGS)
	rm $(OBJECTS)

static_lib: static_lib_no_install
	mv $(STATIC_LIB) $(CGAL_LIB_DESTINATION)

static_lib_no_install: $(OBJECTS)
	$(CGAL_LIB_CREATE)$(STATIC_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_LIB_LDFLAGS) $(STATIC_LIB_ADDITIONNAL_LDFLAGS)
	$(RANLIB) $(STATIC_LIB)
	rm $(OBJECTS)

.PHONY: clean shared_lib shared_lib_no_install status_lib static_lib_no_install

clean:
	rm -f $(STATIC_LIB) $(COMMON_SHARED_LIB) $(OBJECTS)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.SUFFIXES: .cpp

.cpp$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS) $(ADDITIONNAL_CXXFLAGS) -c $<

