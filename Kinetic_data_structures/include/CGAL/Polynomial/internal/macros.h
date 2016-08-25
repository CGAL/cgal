// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_INTERNAL_MACROS_H
#define CGAL_POLYNOMIAL_INTERNAL_MACROS_H

#include <CGAL/Polynomial/internal/config.h>

#ifdef CGAL_POLYNOMIAL_USE_CGAL
/*
  When CGAL is present
*/
#include <CGAL/basic.h>

#define CGAL_POLYNOMIAL_NS CGAL::POLYNOMIAL
#define CGAL_Polynomial_assertion(x) CGAL_assertion(x)
#define CGAL_Polynomial_assertion_code(x) CGAL_assertion_code(x)
#define CGAL_Polynomial_precondition(x) CGAL_precondition(x)
#define CGAL_Polynomial_precondition_code(x) CGAL_precondition_code(x)
#define CGAL_Polynomial_postcondition(x) CGAL_postcondition(x)
#ifdef CGAL_POLYNOMIAL_CHECK_EXPENSIVE
#define CGAL_Polynomial_expensive_precondition(x) CGAL_expensive_precondition(x)
#define CGAL_Polynomial_expensive_assertion(x) CGAL_expensive_assertion(x)
#define CGAL_Polynomial_expensive_postcondition(x) CGAL_expensive_postcondition(x)
#else
#define CGAL_Polynomial_expensive_precondition(x)
#define CGAL_Polynomial_expensive_assertion(x)
#define CGAL_Polynomial_expensive_postcondition(x)
#endif
#define CGAL_Polynomial_exactness_assertion(x) CGAL_exactness_assertion(x)
#define CGAL_Polynomial_exactness_postcondition(x) CGAL_exactness_postcondition(x)
#define CGAL_Polynomial_exactness_precondition(x) CGAL_exactness_precondition(x)

#else
/*
  When no CGAL is present
*/

#define POLYNOMIAL_NS Polynomial

#include <cassert>

#define CGAL_Polynomial_assertion(x) CGAL_assertion(x)
// This does not work
#define CGAL_Polynomial_assertion_code(x) x
#define CGAL_Polynomial_precondition(x) CGAL_assertion(x)
#define CGAL_Polynomial_postcondition(x) CGAL_assertion(x)
#define CGAL_Polynomial_expensive_precondition(x)
#define CGAL_Polynomial_expensive_assertion(x)
#define CGAL_Polynomial_expensive_postcondition(x)
#define CGAL_Polynomial_exactness_postcondition(x)
#define CGAL_Polynomial_exactness_precondition(x)
#endif

#endif
