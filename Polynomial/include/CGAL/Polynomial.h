// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Michael Seel <seel@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@informatik.uni-mainz.de> 
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

/*! \file CGAL/Polynomial.h
 *  \brief Defines class CGAL::Polynomial.
 *  
 *  Polynomials in one variable (or more, by recursion)
 */

#ifndef CGAL_POLYNOMIAL_H
#define CGAL_POLYNOMIAL_H

#include <CGAL/disable_warnings.h>

#include <cstdarg>
#include <cctype>
#include <vector>
#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <boost/mpl/if.hpp>
#include <CGAL/Flattening_iterator.h>

#include <CGAL/Exponent_vector.h>

#include <CGAL/assertions.h>

#ifdef CGAL_USE_LEDA
#  include <LEDA/core/array.h>
#endif // CGAL_USE_LEDA

#include <CGAL/Polynomial/fwd.h>
#include <CGAL/Polynomial/Polynomial_type.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial/Algebraic_structure_traits.h>
#include <CGAL/Polynomial/Real_embeddable_traits.h>
#include <CGAL/Polynomial/Fraction_traits.h>
#include <CGAL/Polynomial/Scalar_factor_traits.h>
#include <CGAL/Polynomial/Modular_traits.h>
#include <CGAL/Polynomial/Coercion_traits.h>
#include <CGAL/Polynomial/Chinese_remainder_traits.h>
#include <CGAL/Polynomial/Get_arithmetic_kernel.h>

// TODO: Are these still includes necessary?
#include <CGAL/Polynomial/polynomial_gcd.h> // used above for Algebraic_structure_traits<Poly...>::Gcd
#include <CGAL/Polynomial/prs_resultant.h>  // for compatibility

#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/polynomial_utils.h>

#include <CGAL/enable_warnings.h>

#endif  // CGAL_POLYNOMIAL_H

// EOF
