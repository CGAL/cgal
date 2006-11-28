# Copyright (c) 2006 Inria Lorraine (France). All rights reserved.
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
#
# Author(s)     : Luis Peñaranda <penarand@loria.fr>

include $(CGAL_MAKEFILE)

EXE_EXT=
LIBMPFI=$(MPFI_LIB)/libmpfi.a
LIBRS=$(RS_LIB)/librs.a
LIBRSE=$(RS_LIB)/librsexport.a
LIBCGALGBRS=$(CGAL_LIB_DESTINATION)/libCGALGbrs.a
LIBS_FLAGS=-lmpfi -lrs -lrsexport -lCGALGbrs
MYFLAGS=-g #-DCGAL_RS_DEBUG -DCGAL_RS_VERB=0

OBJS= \
      Gbrs_polynomial_1$(EXE_EXT) \
      Gbrs_polynomial_2$(EXE_EXT) \
      Gbrs_algebraic_kernel_1$(EXE_EXT)

outside: OUTSIDEPARS=-I$(PWD)/../../include -I$(RS_INCL) -L$(RS_LIB)
outside: $(OBJS)

all: $(OBJS)

clean:
	rm -f $(OBJS)

.SUFFIXES: .C $(EXE_EXT) $(SUFFIXES)

.C$(EXE_EXT):
	$(CGAL_CXX) $(MYFLAGS) $(OUTSIDEPARS) $(CGAL_CXXFLAGS) $(CGAL_LIBPATHFLAGS) $< $(CGAL_LDFLAGS) $(LIBS_FLAGS) $(LIBMPFI) $(LIBRSE) $(LIBRS) $(LIBCGALGBRS) -o $@

