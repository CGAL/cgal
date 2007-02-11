# Copyright (c) 2006 Tel-Aviv University (Israel)
# Copyright (c) 2007 INRIA Sophia-Antipolis (France)
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
# Author(s): Efi Fogel         <efif@post.tau.ac.il>
#            Laurent Rineau

# defaut rule
defaut: static_lib

#---------------------------------------------------------------------#
#                    object files
#---------------------------------------------------------------------#

OBJECTS = \
Benchmark$(OBJ_EXT) \
benchmark_parser$(OBJ_EXT) \
benchmark_lexer$(OBJ_EXT) \
Benchmark_visitor$(OBJ_EXT) \
Option_parser$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    libCGAL settings
#---------------------------------------------------------------------#

SOVERSION=
SOMAJOR=

SHARED_LIB_CREATE = $(CGAL_BENCHMARK_SHARED_LIB_CREATE)
SHARED_LIB        = $(CGAL_BENCHMARK_SHARED_LIB)
STATIC_LIB        = $(CGAL_BENCHMARK_LIB)

#---------------------------------------------------------------------#
#                    common rules for all libraries
#---------------------------------------------------------------------#
include ../CGAL/common.mk

