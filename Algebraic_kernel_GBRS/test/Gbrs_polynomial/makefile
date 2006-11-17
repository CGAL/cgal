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

include ${CGAL_MAKEFILE}

MYFLAGS=-g #-DCGAL_RS_DEBUG -DCGAL_RS_VERB=0

OBJS=Gbrs_polynomial_1 Gbrs_polynomial_2 Gbrs_algebraic_kernel_1

all: ${OBJS}

# TODO: make general rules

Gbrs_polynomial_1: Gbrs_polynomial_1.C
	${CGAL_CXX} ${MYFLAGS} ${CGAL_CXXFLAGS} ${CGAL_LIBPATHFLAGS} $< ${CGAL_LDFLAGS} -o $@

Gbrs_polynomial_2: Gbrs_polynomial_2.C
	${CGAL_CXX} ${MYFLAGS} ${CGAL_CXXFLAGS} ${CGAL_LIBPATHFLAGS} $< ${CGAL_LDFLAGS} -o $@

Gbrs_algebraic_kernel_1: Gbrs_algebraic_kernel_1.C
	${CGAL_CXX} ${MYFLAGS} ${CGAL_CXXFLAGS} ${CGAL_LIBPATHFLAGS} $< ${CGAL_LDFLAGS} -o $@

clean:
	rm -f ${OBJS}

