# Copyright (c) 2004,2006-2007 INRIA Sophia-Antipolis (France).
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



### inclusion of KDS sub-makefile

# This rule to be able to compile Qt_widget without KDS.
# For Solaris make, need to be *before* the inclusion.
kds.makefile:
	touch kds.makefile

include kds.makefile

#---------------------------------------------------------------------#
#                    object files
#---------------------------------------------------------------------#

OBJECTS = \
Qt_widget$(OBJ_EXT) \
Qt_widget_layer$(OBJ_EXT) \
Qt_widget_standard_toolbar$(OBJ_EXT) \
Qt_help_window$(OBJ_EXT) \
Qt_widget_history$(OBJ_EXT) \
Qt_widget_xpm_icons$(OBJ_EXT) \
Qt_widget_OpenGL$(OBJ_EXT) \
$(CGAL_KDS_OBJECTS)


### inclusion of KDS dependencies sub-makefile

# This rule to be able to compile Qt_widget without KDS:
# For Solaris make, need to be *before* the inclusion.
kds_deps.makefile:
	touch kds_deps.makefile

include kds_deps.makefile


#---------------------------------------------------------------------#
#                    specific rules for libCGALQt
#---------------------------------------------------------------------#

Qt_widget$(OBJ_EXT):
	$(QT_MOC) '$(CGAL_INCL_DIR)/CGAL/IO/Qt_widget.h' -o Qt_widget.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget.cpp

Qt_widget_layer$(OBJ_EXT):
	$(QT_MOC) '$(CGAL_INCL_DIR)/CGAL/IO/Qt_widget_layer.h' -o Qt_widget_layer.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_layer.cpp

Qt_widget_standard_toolbar$(OBJ_EXT):
	$(QT_MOC) '$(CGAL_INCL_DIR)/CGAL/IO/Qt_widget_standard_toolbar.h' -o Qt_widget_standard_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_standard_toolbar.cpp

Qt_help_window$(OBJ_EXT):
	$(QT_MOC) '$(CGAL_INCL_DIR)/CGAL/IO/Qt_help_window.h' -o Qt_help_window.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_help_window.cpp

Qt_widget_history$(OBJ_EXT):
	$(QT_MOC) '$(CGAL_INCL_DIR)/CGAL/IO/Qt_widget_history.h' -o Qt_widget_history.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_history.cpp

Qt_widget_xpm_icons$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_xpm_icons.cpp

Qt_widget_OpenGL$(OBJ_EXT):
	$(QT_MOC) '$(CGAL_INCL_DIR)/CGAL/IO/Qt_widget_OpenGL.h' -o Qt_widget_OpenGL.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_OpenGL.cpp

#---------------------------------------------------------------------#
#                    libCGALQt settings
#---------------------------------------------------------------------#

SOVERSION=2.0.0
SOMAJOR=2

SHARED_LIB_CREATE = $(CGAL_QT_SHARED_LIB_CREATE)
SHARED_LIB        = $(CGAL_QT_SHARED_LIB)
STATIC_LIB        = $(CGAL_QT_LIB)

#---------------------------------------------------------------------#
#                    common rules for all libraries
#---------------------------------------------------------------------#
include ../CGAL/common.mk

