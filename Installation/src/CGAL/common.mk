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
#                    target entries
#---------------------------------------------------------------------#

shared_lib: shared_lib_no_install
	mv $(CGAL_SHARED_LIB_WITH_SOVERSION) $(CGAL_LIB_DESTINATION)
	if [ "$(CGAL_SHARED_LIB_WITH_SOVERSION)" != "$(SHARED_LIB)" ]; then \
	  rm -f $(CGAL_LIB_DESTINATION)/"$(SHARED_LIB)"; \
	  ln -s "$(CGAL_SHARED_LIB_WITH_SOVERSION)" $(CGAL_LIB_DESTINATION)/"$(SHARED_LIB)"; \
	fi
	if [ "$(CGAL_SHARED_LIB_WITH_SOVERSION)" != "$(CGAL_SHARED_LIB_WITH_SOMAJOR)" ]; then \
	  rm -f $(CGAL_LIB_DESTINATION)/"$(CGAL_SHARED_LIB_WITH_SOMAJOR)"; \
	  ln -s "$(CGAL_SHARED_LIB_WITH_SOVERSION)" $(CGAL_LIB_DESTINATION)/"$(CGAL_SHARED_LIB_WITH_SOMAJOR)"; \
	fi

shared_lib_no_install: $(OBJECTS)
	$(CGAL_SHARED_LIB_CREATE)$(CGAL_SHARED_LIB_WITH_SOVERSION) $(CGAL_SHARED_LIB_SONAME) \
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
	rm -f $(STATIC_LIB) $(CGAL_SHARED_LIB_WITH_SOVERSION) $(OBJECTS)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

%$(OBJ_EXT) : %.cpp
	$(CGAL_CXX)  $(CXXFLAGS) $(ADDITIONNAL_CXXFLAGS) -c $<

FLEX ?=flex
BISON ?=bison

%.cpp : %.l
	$(FLEX) -8 -o$@ $<

# Apparently, old versions of bison (e.g., 1.28) name the generated definion
# header file <base>.cpp.h. The file must be renamed to <base>.hpp
BISON_EXISTS_CMD =which $(BISON)
BISON_EXISTS =$(shell $(BISON_EXISTS_CMD))
ifneq ($(strip $(BISON_EXISTS)),)
  BISON_VERSION_CMD =expr match "`$(BISON) --version`" '.*\([1-9]\.[0-9]*\)'
  BISON_VERSION =$(shell $(BISON_VERSION_CMD))
  OLD_BISON_VERSION_CMD =expr "$(BISON_VERSION)" \<= 1.28
  OLD_BISON_VERSION =$(shell $(OLD_BISON_VERSION_CMD))
endif

%.cpp %.hpp : %.y
	$(BISON) -d $< -o $*.cpp
ifeq ($(OLD_BISON_VERSION), 1)
	mv $*.cpp.h $*.hpp
endif
