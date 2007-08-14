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

default-qtwidget: $(CGAL_SHARED_LIBNAME)

#---------------------------------------------------------------------#
#                    KDS sub-makefile
#---------------------------------------------------------------------#

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

#---------------------------------------------------------------------#
#                    libCGALQt settings
#---------------------------------------------------------------------#

SOVERSION=2.0.1
SOMAJOR=2

LIB_NAME_BASE = CGALQt

ADDITIONAL_CXXFLAGS = $(ALL_QT_CXXFLAGS) $(ALL_OPENGL_CXXFLAGS) -I.
SHARED_LIB_ADDITIONAL_LDFLAGS = \
            $(QT_LIBPATHFLAGS) $(QT_LDFLAGS) \
            $(OPENGL_LIBPATHFLAGS) $(OPENGL_LDFLAGS) \
            $(X11_LIBPATHFLAGS) $(X11_LDFLAFS) \
            $(CGAL_LIBPATHFLAGS) $(CGAL_CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    common rules for all libraries
#---------------------------------------------------------------------#
CGAL_COMMON_MK = ../CGAL/common.mk
include $(CGAL_COMMON_MK)

#---------------------------------------------------------------------#
#                    specific rules for libCGALQt
#---------------------------------------------------------------------#

%.moc : %.h
	$(QT_MOC) $< -o $@

vpath %.h ../../include/CGAL/IO
vpath %.h ../../../Nef_S2/include/CGAL/IO/
vpath Kinetic_% ../../../Kinetic_data_structures/src/CGALQt
vpath Qt_widget_OpenGL.cpp ../../../Nef_S2/src/CGALQt/

clean_temp_files::
	rm -f *.moc

#---------------------------------------------------------------------#
#                    moc rules
#---------------------------------------------------------------------#

Qt_widget$(OBJ_EXT) : Qt_widget.moc

Qt_widget_layer$(OBJ_EXT) : Qt_widget_layer.moc

Qt_widget_standard_toolbar$(OBJ_EXT) : Qt_widget_standard_toolbar.moc

Qt_help_window$(OBJ_EXT) : Qt_help_window.moc

Qt_widget_history$(OBJ_EXT) : Qt_widget_history.moc

Qt_widget_OpenGL$(OBJ_EXT): Qt_widget_OpenGL.moc

#---------------------------------------------------------------------#
#                    KDS sub-makefile
#---------------------------------------------------------------------#

# This rule to be able to compile Qt_widget without KDS.
# For Solaris make, need to be *before* the inclusion.
kds_deps.makefile:
	touch kds_deps.makefile

include kds_deps.makefile
