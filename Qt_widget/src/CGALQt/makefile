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

vpath %.cpp ../../../Kinetic_data_structures/src/CGALQt

# This rule to be able to compile Qt_widget without KDS.
# For Solaris make, need to be *before* the inclusion.
kds.makefile:
	touch kds.makefile

ifneq ($(wildcard kds.makefile),)
  include kds.makefile
else
  ifneq ($(wildcard ../../../Kinetic_data_structures/src/CGALQt/kds.makefile),)
    include ../../../Kinetic_data_structures/src/CGALQt/kds.makefile
  endif
endif

#---------------------------------------------------------------------#
#                    object files
#---------------------------------------------------------------------#

# We now compile all source files at once since it is faster
# When adding a new file, please update app_files.cpp
OBJECTS = all_files$(OBJ_EXT)

# Below variable is not used now (but may come back so let's keep it)
OBJECTS_OLD = \
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

SOVERSION=2.0.0
SOMAJOR=2

LIB_NAME_BASE = CGALQt

ADDITIONAL_CXXFLAGS = $(ALL_QT_CXXFLAGS) $(ALL_OPENGL_CXXFLAGS) \
			-I. \
			-I ../../../Kinetic_data_structures/src/CGALQt \
			-I ../../../Nef_S2/src/CGALQt/
SHARED_LIB_ADDITIONAL_LDFLAGS = \
            $(QT_LIBPATHFLAGS) $(QT_LDFLAGS) \
            $(OPENGL_LIBPATHFLAGS) $(OPENGL_LDFLAGS) \
            $(X11_LIBPATHFLAGS) $(X11_LDFLAFS) \
            $(CGAL_LIBPATHFLAGS) $(CGAL_CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    common rules for all libraries
#---------------------------------------------------------------------#
CGAL_COMMON_MK = ../CGAL/common.mk
ifeq ($(willcard $(CGAL_COMMON_MK)),)
  CGAL_COMMON_MK = ../../../Installation/src/CGAL/common.mk
endif
include $(CGAL_COMMON_MK)

#---------------------------------------------------------------------#
#                    specific rules for libCGALQt
#---------------------------------------------------------------------#

%.moc : %.h
	$(QT_MOC) $< -o $@

vpath %.h ../../include/CGAL/IO ../../../Nef_S2/include/CGAL/IO/ ../../../Kinetic_data_structures/include/CGAL/Kinetic/IO/internal ../../include/CGAL/Kinetic/IO/internal
vpath Kinetic_% ../../../Kinetic_data_structures/src/CGALQt
vpath Qt_widget_OpenGL.cpp ../../../Nef_S2/src/CGALQt/

clean_temp_files::
	rm -f *.moc

.PHONY: regenerate_all_files.cpp

regenerate_all_files:
	echo $(OBJECTS_OLD)
	rm -f all_files.cpp
	touch all_files.cpp
	for f in $(foreach file,$(OBJECTS_OLD),$(basename $(file))); do \
	  echo "#include \"$$f.cpp\"" >> all_files.cpp; \
	done

#---------------------------------------------------------------------#
#                    moc rules
#---------------------------------------------------------------------#

Qt_widget$(OBJ_EXT) : Qt_widget.moc

Qt_widget_layer$(OBJ_EXT) : Qt_widget_layer.moc

Qt_widget_standard_toolbar$(OBJ_EXT) : Qt_widget_standard_toolbar.moc

Qt_help_window$(OBJ_EXT) : Qt_help_window.moc

Qt_widget_history$(OBJ_EXT) : Qt_widget_history.moc

Qt_widget_OpenGL$(OBJ_EXT): Qt_widget_OpenGL.moc

all_files$(OBJ_EXT) : \
	Qt_widget.moc \
	Qt_widget_layer.moc \
	Qt_widget_standard_toolbar.moc \
	Qt_help_window.moc \
	Qt_widget_history.moc \
	Qt_widget_OpenGL.moc \
	Kinetic_Qt_widget_2_core.moc \
	Kinetic_Qt_core.moc \
	Kinetic_Qt_timer.moc \
	Kinetic_Qt_window_2.moc \
	$(wildcard Kinetic_*.xpm)

#---------------------------------------------------------------------#
#                    KDS sub-makefile
#---------------------------------------------------------------------#

# This rule to be able to compile Qt_widget without KDS.
# For Solaris make, need to be *before* the inclusion.
kds_deps.makefile:
	touch kds_deps.makefile

ifneq ($(wildcard kds_deps.makefile),)
  include kds_deps.makefile
else
  ifneq ($(wildcard ../../../Kinetic_data_structures/src/CGALQt/kds_deps.makefile),)
    include ../../../Kinetic_data_structures/src/CGALQt/kds_deps.makefile
  endif
endif
