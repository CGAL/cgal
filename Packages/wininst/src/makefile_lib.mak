# Copyright (c) 1999,2001 The CGAL Consortium
#
# Package: wininst
#
# Used for Windows-specific installation
#
#
# This software and related documentation is part of the
# Computational Geometry Algorithms Library (CGAL).
#
# This is the nmake makefile for compiling the CGAL object library 
# CGAL.lib using MSVC compiler
#
#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#

!INCLUDE $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = $(CGAL_LIB_CXXFLAGS)

#---------------------------------------------------------------------#
# Object files to create
#---------------------------------------------------------------------#

CGAL_OBJECTS = \
Bbox_2$(OBJ_EXT) \
Bbox_2_intersections$(OBJ_EXT) \
Bbox_3_intersections$(OBJ_EXT) \
Color$(OBJ_EXT) \
File_header_OFF$(OBJ_EXT) \
File_header_extended_OFF$(OBJ_EXT) \
File_scanner_OFF$(OBJ_EXT) \
File_writer_OFF$(OBJ_EXT) \
File_writer_VRML_2$(OBJ_EXT) \
File_writer_inventor$(OBJ_EXT) \
File_writer_wavefront$(OBJ_EXT) \
Geomview_stream$(OBJ_EXT) \
Interval_arithmetic$(OBJ_EXT) \
MP_Float$(OBJ_EXT) \
Origin$(OBJ_EXT) \
Random$(OBJ_EXT) \
Triangulation_3$(OBJ_EXT) \
aff_transformation_tags$(OBJ_EXT) \
assertions$(OBJ_EXT) \
io$(OBJ_EXT) \
optimisation_basic$(OBJ_EXT)

#---------------------------------------------------------------------#
# Files to store in the library
#---------------------------------------------------------------------#

CGAL_OBJECTS_LIBPARAM = \
$(CGAL_OBJ_PREFIX)Bbox_2$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)Bbox_2_intersections$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)Bbox_3_intersections$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)Color$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)File_header_OFF$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)File_header_extended_OFF$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)File_scanner_OFF$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)File_writer_OFF$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)File_writer_VRML_2$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)File_writer_inventor$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)File_writer_wavefront$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)Geomview_stream$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)Interval_arithmetic$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)MP_Float$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)Origin$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)Random$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)Triangulation_3$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)aff_transformation_tags$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)assertions$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)io$(OBJ_EXT) \
$(CGAL_OBJ_PREFIX)optimisation_basic$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

lib: lib_no_install
	$(MV) $(CGAL_LIB) $(CGAL_LIB_DESTINATION)\$(CGAL_LIB)

lib_no_install: $(CGAL_OBJECTS) $(CGAL_EXTRA_OBJECTS)
	$(CGAL_LIB_CREATE)$(CGAL_LIB) \
	$(CGAL_OBJECTS_LIBPARAM) \
	$(CGAL_OBJ_PREFIX)$(CGAL_EXTRA_OBJECTS) $(CGAL_LIB_LDFLAGS)
	$(RM) *.obj

#	$(RM) $(CGAL_OBJECTS) $(CGAL_EXTRA_OBJECTS)
# 		does not work on Win95



workaround_4_ms$(OBJ_EXT):
	$(CP) Interval_arithmetic\workaround_4_ms.obj workaround_4_ms$(OBJ_EXT)

clean:
	$(RM_FORCE) $(CGAL_LIB) $(CGAL_OBJECTS) $(CGAL_EXTRA_OBJECTS)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS)  -c $<

