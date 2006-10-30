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

#ifndef CGAL_RS_VERB
#ifdef CGAL_RS_DEBUG
#define CGAL_RS_VERB 1
#else
#define CGAL_RS_VERB 0
#endif
#endif

// the default precision of RS to calculate a root (precision is 2^n)
#ifndef CGAL_RS_DEF_PREC
#define CGAL_RS_DEF_PREC 5
#endif

// the minimum, used when calculating a sign
#ifndef CGAL_RS_MIN_PREC
#define CGAL_RS_MIN_PREC 5
#endif

// when refining a calculation, increase by this factor
#ifndef CGAL_RS_PREC_FACTOR
#define CGAL_RS_PREC_FACTOR 2
#endif

// after reaching this precision, give up
#ifndef CGAL_RS_MAX_PREC
#define CGAL_RS_MAX_PREC 80
#endif

#include <mpfi.h>
#include <CGAL/Gbrs_polynomial_1.h>
#include <CGAL/Gbrs_algebraic_1.h>

CGAL_BEGIN_NAMESPACE

// initialize the RS solver, returns 0 if everything was OK
int init_solver ();

// solve given the precision, returns de number of roots
int solve_1 (mpfi_ptr *&, const Rational_polynomial_1 &,
		unsigned int = CGAL_RS_DEF_PREC);

// evaluate a polynomial at a given algebraic number
Sign sign_1 (const Rational_polynomial_1 &, const Algebraic_1 &);

// compare two algebraic numbers
Comparison_result compare_1 (Algebraic_1 &, Algebraic_1 &);

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Gbrs_solve_1.C>
#endif	// CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif	// CGAL_GBRS_SOLVE_1_H
