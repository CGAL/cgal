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

# CXXFLAGS can be passed as command line argument
CXXFLAGS = $(CGAL_SHARED_LIB_CXXFLAGS)


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

.PHONY: default

default: $(CGAL_SHARED_LIBNAME)

.PHONY: shared_lib_no_install shared_lib

LIBRARIES_EXTRA_NAMES :=
ifneq "$(CGAL_SHARED_LIBNAME_WITH_SOVERSION)" "$(CGAL_SHARED_LIBNAME)"
  LIBRARIES_EXTRA_NAMES += $(CGAL_SHARED_LIBNAME)
endif
ifneq "$(CGAL_SHARED_LIBNAME_WITH_SOMAJOR)" "$(CGAL_SHARED_LIBNAME)"
  ifneq "$(CGAL_SHARED_LIBNAME_WITH_SOMAJOR)" "$(CGAL_SHARED_LIBNAME_WITH_SOVERSION)"
    LIBRARIES_EXTRA_NAMES += $(CGAL_SHARED_LIBNAME_WITH_SOMAJOR)
  endif
endif

ifneq "$(LIBRARIES_EXTRA_NAMES)" ""
###### If CGAL_SHARED_LIBNAME_WITH_SOMAJOR is different from CGAL_SHARED_LIBNAME.

shared_lib_no_install: $(CGAL_SHARED_LIBNAME_WITH_SOVERSION)
	$(MAKE) clean_temp_files

shared_lib: shared_lib_no_install
	mv $(CGAL_SHARED_LIBNAME_WITH_SOVERSION) $(CGAL_LIB_DESTINATION)
	for symlink in $(LIBRARIES_EXTRA_NAMES); do \
	  mv "$$symlink" $(CGAL_LIB_DESTINATION); \
	done

$(CGAL_SHARED_LIBNAME_WITH_SOVERSION) $(LIBRARIES_EXTRA_NAMES): $(OBJECTS)
	$(CGAL_SHARED_LIB_CREATE)$(CGAL_SHARED_LIBNAME_WITH_SOVERSION) $(CGAL_SHARED_LIB_SONAME) \
	  $(OBJECTS:%=$(CGAL_OBJ_PREFIX)%) \
	  $(CGAL_SHARED_LIB_LDFLAGS) $(SHARED_LIB_ADDITIONAL_LDFLAGS)
	for symlink in $(LIBRARIES_EXTRA_NAMES); do \
	  rm -f "$$symlink"; \
	  ln -s "$(CGAL_SHARED_LIBNAME_WITH_SOVERSION)" "$$symlink"; \
	done
else
###### If CGAL_SHARED_LIBNAME_WITH_SOVERSION is the same as CGAL_SHARED_LIBNAME.
shared_lib_no_install: $(CGAL_SHARED_LIBNAME)
	$(MAKE) clean_temp_files

shared_lib: shared_lib_no_install
	mv $(CGAL_SHARED_LIBNAME) $(CGAL_LIB_DESTINATION)

$(CGAL_SHARED_LIBNAME): $(OBJECTS)
	$(CGAL_SHARED_LIB_CREATE)$(CGAL_SHARED_LIBNAME) \
	  $(OBJECTS:%=$(CGAL_OBJ_PREFIX)%) \
	  $(CGAL_SHARED_LIB_LDFLAGS) $(SHARED_LIB_ADDITIONAL_LDFLAGS)
endif

.PHONY: static_lib_no_install static_lib

static_lib_no_install: $(CGAL_STATIC_LIBNAME)
	$(MAKE) clean_temp_files

static_lib: static_lib_no_install
	mv $(CGAL_STATIC_LIBNAME) $(CGAL_LIB_DESTINATION)

$(CGAL_STATIC_LIBNAME): $(OBJECTS)
	$(CGAL_STATIC_LIB_CREATE)$(CGAL_STATIC_LIBNAME) \
	  $(OBJECTS:%=$(CGAL_OBJ_PREFIX)%) \
	  $(CGAL_LIB_LDFLAGS) $(STATIC_LIB_ADDITIONAL_LDFLAGS)
	$(RANLIB) $(CGAL_STATIC_LIBNAME)

clean_temp_files::
	rm -f $(OBJECTS)

clean:: clean_temp_files
	rm -f $(CGAL_STATIC_LIBNAME)
	rm -f $(CGAL_SHARED_LIBNAME) $(LIBRARIES_EXTRA_NAMES)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

%$(OBJ_EXT): %.cpp
	$(CGAL_CXX)  $(CXXFLAGS) $(ADDITIONAL_CXXFLAGS) -c $<

FLEX ?=flex
BISON ?=bison

%.cpp: %.l
	$(FLEX) -8 -o$@ $<

%.cpp %.hpp: %.y
	$(BISON) -d $< -o $*.cpp
# Apparently, old versions of bison (e.g., 1.28) name the generated definion
# header file <base>.cpp.h. The file must be renamed to <base>.hpp
	@if [ ! -e $*.hpp -a -e $*.cpp.h ]; then \
	  echo Moving $*.hpp to $*.cpp.h; \
	  mv $*.cpp.h $*.hpp; \
	fi
