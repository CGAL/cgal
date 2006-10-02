// Copyright (c) 2006 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Luis Peñaranda <penarand@loria.fr>

#ifndef CGAL_GBRS_SOLVE_1_H
#define CGAL_GBRS_SOLVE_1_H

#include <mpfi.h>
#include <CGAL/Gbrs_polynomial_1.h>

CGAL_BEGIN_NAMESPACE
// both functions return the number of roots, -1 if it is infinite and -2 if
// there were an error

// solve given the precision
int solve_1 (mpfi_t *&, const Rational_polynomial_1 &, unsigned int);

// solve with the default precision
int solve_1 (mpfi_t *&, const Rational_polynomial_1 &);

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Gbrs_solve_1.cpp>
#endif	// CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif	// CGAL_GBRS_SOLVE_1_H
