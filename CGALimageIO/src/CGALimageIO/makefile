# Copyright (c) 2007 INRIA Sophia-Antipolis (France).
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

OBJECTS = \
  Image_3$(OBJ_EXT) \
  all_files$(OBJ_EXT)
#   ImageIO$(OBJ_EXT) \
#   analyze$(OBJ_EXT) \
#   bmp$(OBJ_EXT) \
#   bmpendian$(OBJ_EXT) \
#   bmpread$(OBJ_EXT) \
#   convert$(OBJ_EXT) \
#   gif$(OBJ_EXT) \
#   gis$(OBJ_EXT) \
#   inr$(OBJ_EXT) \
#   iris$(OBJ_EXT) \
#   pnm$(OBJ_EXT) \
#   recbuffer$(OBJ_EXT) \
#   recline$(OBJ_EXT) \
#   reech4x4$(OBJ_EXT) \
#   fgetns$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    libCGALimageIO settings
#---------------------------------------------------------------------#

SOVERSION=3.0.0
SOMAJOR=3

LIB_NAME_BASE = CGALimageIO

ADDITIONAL_CXXFLAGS = -I../../include $(OTHER_CXXFLAGS) $(ALL_ZLIB_CXXFLAGS) $(ALL_OPENGL_CXXFLAGS) 
SHARED_LIB_ADDITIONAL_LDFLAGS = $(ZLIB_LIBPATHFLAGS) $(ZLIB_LDFLAGS) $(CGAL_LIBPATHFLAGS)  $(CGAL_CGAL_LDFLAGS) $(OPENGL_LIBPATHFLAGS) $(OPENGL_LDFLAGS)

#---------------------------------------------------------------------#
#                    common rules for all libraries
#---------------------------------------------------------------------#
include ../CGAL/common.mk


