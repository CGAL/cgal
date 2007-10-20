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
Option_parser$(OBJ_EXT) \
# Polynomial_visitor$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    libCGAL settings
#---------------------------------------------------------------------#

SOVERSION=
SOMAJOR=
LIB_NAME_BASE = CGALBenchmark

# The following is needed for some linkers (e.g., Darwin)
SHARED_LIB_ADDITIONAL_LDFLAGS = \
            $(BOOST_PROGRAM_OPTIONS_LIBPATHFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS)

#---------------------------------------------------------------------#
#                    common rules for all libraries
#---------------------------------------------------------------------#
include ../CGAL/common.mk

#---------------------------------------------------------------------#
#                    additional rules
#---------------------------------------------------------------------#

clean_temp_files::
	rm -f benchmark_parser.hpp benchmark_lexer.cpp benchmark_parser.cpp

benchmark_lexer$(OBJ_EXT): benchmark_parser.hpp
