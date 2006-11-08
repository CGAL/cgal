# Copyright (c) 2004,2006  INRIA Sophia-Antipolis (France).
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
# Author(s)     : various

# This is the makefile for compiling the CGAL Qt support library libCGALQt.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = SHOULD_BE_SET_BY_COMMAND_LINE
include $(CGAL_MAKEFILE) 

# defaut rule
all: lib

### inclusion of KDS sub-makefile

# This rule to be able to compile Qt_widget without KDS.
# For Solaris make, need to be *before* the inclusion.
kds.makefile:
	touch kds.makefile

include kds.makefile

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

# CXXFLAGS must be provided on the make command-line.
# CXXFLAGS = -I"../../include" $(CGAL_WINLIB_CXXFLAGS) -DCGAL_USE_QT

#---------------------------------------------------------------------#
# Object files
#---------------------------------------------------------------------#

CGAL_OBJECTS = \
Qt_widget$(OBJ_EXT) \
Qt_widget_layer$(OBJ_EXT) \
Qt_widget_standard_toolbar$(OBJ_EXT) \
Qt_help_window$(OBJ_EXT) \
Qt_widget_history$(OBJ_EXT) \
Qt_widget_xpm_icons$(OBJ_EXT) \
Qt_widget_OpenGL$(OBJ_EXT) \
$(CGAL_KDS_OBJECTS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

shared_lib: $(CGAL_OBJECTS)
	$(CGAL_SHARED_LIB_CREATE)$(CGAL_QT_SHARED_LIB) -Wl,-soname,$(CGAL_QT_SHARED_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_SHARED_LIB_LDFLAGS) -lGL
	mv $(CGAL_QT_SHARED_LIB) '$(CGAL_LIB_DESTINATION)'
	rm $(CGAL_OBJECTS)

static_lib: lib

lib: lib_no_install
	mv $(CGAL_QT_LIB) $(CGAL_LIB_DESTINATION)

lib_no_install: $(CGAL_OBJECTS) 
	$(CGAL_LIB_CREATE)$(CGAL_QT_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_LIB_LDFLAGS)
	$(RANLIB) $(CGAL_QT_LIB)
	rm $(CGAL_OBJECTS) 

clean:
	rm -f $(CGAL_QT_LIB) $(CGAL_QT_SHARED_LIB) $(CGAL_OBJECTS) *.moc

### inclusion of KDS dependencies sub-makefile

# This rule to be able to compile Qt_widget without KDS:
# For Solaris make, need to be *before* the inclusion.
kds_deps.makefile:
	touch kds_deps.makefile

include kds_deps.makefile

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

Qt_widget$(OBJ_EXT):
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/IO/Qt_widget.h -o Qt_widget.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget.cpp

Qt_widget_layer$(OBJ_EXT):
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/IO/Qt_widget_layer.h -o Qt_widget_layer.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_layer.cpp

Qt_widget_standard_toolbar$(OBJ_EXT):
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/IO/Qt_widget_standard_toolbar.h -o Qt_widget_standard_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_standard_toolbar.cpp

Qt_help_window$(OBJ_EXT):
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/IO/Qt_help_window.h -o Qt_help_window.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_help_window.cpp

Qt_widget_history$(OBJ_EXT):
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/IO/Qt_widget_history.h -o Qt_widget_history.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_history.cpp

Qt_widget_xpm_icons$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_xpm_icons.cpp

Qt_widget_OpenGL$(OBJ_EXT):
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/IO/Qt_widget_OpenGL.h -o Qt_widget_OpenGL.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_OpenGL.cpp

.SUFFIXES: .cpp $(SUFFIXES)

.cpp$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS) -c $<

