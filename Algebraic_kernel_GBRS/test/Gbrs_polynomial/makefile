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
#LIBRS=$(RS_LIB)/librs.a
LIBRS=
LIBRSE=$(RS_LIB)/librsexport.a
LIBCGALGBRS=$(CGAL_LIB_DESTINATION)/libCGALGbrs.a
#LIBS_FLAGS=-lmpfi -lCGALGbrs
#LIBS_FLAGS=-lmpfi -lrsexport -lCGALGbrs
LIBS_FLAGS=-lmpfi -lrsexport -lCGALGbrs
MYFLAGS=-O3 -D__TEST_ARR #-g -O0 #-DUSE_FGB #-DCGAL_RS_DEBUG -DCGAL_RS_VERB=0


DISTRIB=${HOME}/Distrib/IBM_INTEL_LINUX.4.2.1
#DISTRIB=/usr/local/science/Distrib/IBM_INTEL_LINUX.4.2.1
#MPFI_INCL=${DISTRIB}/include
RS_INCL=${DISTRIB}/include
#MPFI_LIB=${DISTRIB}/lib
RS_LIB=${DISTRIB}/lib
LIB_STAT_PATH=$(RS_LIB)
#FGBFLAGS=-DUSE_FGB -L${LIB_STAT_PATH} ${LIB_STAT_PATH}/libfgb2rs.a ${LIB_STAT_PATH}/libfgbexp.a ${LIB_STAT_PATH}/libfgb.a ${LIB_STAT_PATH}/libgbvgf.a ${LIB_STAT_PATH}/libminpolyvgf.a
#FGBFLAGS=-L${LIB_STAT_PATH} ${LIB_STAT_PATH}/libfgb2rs.a ${LIB_STAT_PATH}/libfgbexp.a ${LIB_STAT_PATH}/libfgb.a ${LIB_STAT_PATH}/libgbvgf.a ${LIB_STAT_PATH}/libminpolyvgf.a
FGBFLAGS=-L${LIB_STAT_PATH} 


OBJS= \
      Gbrs_polynomial_1$(EXE_EXT) \
      Gbrs_algebraic_kernel_1$(EXE_EXT) \
      Arr_poly_traits_1$(EXE_EXT) \
      Arrangement$(EXE_EXT) \
      Core$(EXE_EXT) \
      Parser_polys$(EXE_EXT) \
      Gbrs_polynomial_2$(EXE_EXT) #\
      Gbrs_solve_2$(EXE_EXT) \

outside: OUTSIDEPARS=-I$(PWD)/../../include -I$(RS_INCL) -L$(RS_LIB)
outside: $(OBJS)

all: $(OBJS)

clean:
	rm -f $(OBJS)

.SUFFIXES: .C $(EXE_EXT) $(SUFFIXES)

.C$(EXE_EXT):
	$(CGAL_CXX) $(MYFLAGS) $(OUTSIDEPARS) $(CGAL_CXXFLAGS) $(CGAL_LIBPATHFLAGS) $< $(CGAL_LDFLAGS) $(LIBMPFI) $(LIBRSE) $(LIBRS) $(LIBCGALGBRS) $(LIBS_FLAGS) $(FGBFLAGS) -o $@

