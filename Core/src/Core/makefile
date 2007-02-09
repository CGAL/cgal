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
CoreAux$(OBJ_EXT) \
extLong$(OBJ_EXT) \
CoreDefs$(OBJ_EXT) \
CoreIO$(OBJ_EXT) \
BigFloat$(OBJ_EXT) \
Real$(OBJ_EXT) \
Expr$(OBJ_EXT) \
GmpIO$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    libCore settings
#---------------------------------------------------------------------#

SOVERSION=2.0.0
SOMAJOR=2

SHARED_LIB_CREATE = $(CGAL_CORE_SHARED_LIB_CREATE)
SHARED_LIB        = $(CGAL_CORE_SHARED_LIB)
STATIC_LIB        = $(CGAL_CORE_LIB)
SHARED_LIB_ADDITIONNAL_LDFLAGS = $(GMP_LIBPATHFLAGS) $(MPFR_LIBPATHFLAGS) $(GMP_LDFLAGS)

# The shared_lib case is still buggy:
# - Currently MPFR is linked with it (through GMP_LDFLAGS), although it is
#   not needed, but anyway next Core release will need it...
# - When GMP or MPFR is a static library on Darwin, it fails because -fno-common is required...

#---------------------------------------------------------------------#
#                    common rules for all libraries
#---------------------------------------------------------------------#
include ../CGAL/common.mk
