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


#---------------------------------------------------------------------#
#                    object files
#---------------------------------------------------------------------#

# We now compile all source files at once since it is faster
# When adding a new file, please update app_files.cpp
OBJECTS = all_files$(OBJ_EXT)

# Below variable is not used now (but may come back so let's keep it)
OBJECTS_OLD = \
assertions$(OBJ_EXT) \
Bbox_2_intersections$(OBJ_EXT) \
Bbox_3_intersections$(OBJ_EXT) \
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
#                    libCGAL settings
#---------------------------------------------------------------------#

SOVERSION=2.0.0
SOMAJOR=2

LIB_NAME_BASE = CGAL

#---------------------------------------------------------------------#
#                    common rules for all libraries
#---------------------------------------------------------------------#
include ../CGAL/common.mk

