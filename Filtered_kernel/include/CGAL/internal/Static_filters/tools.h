// Copyright (c) 2001,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_INTERNAL_STATIC_FILTERS_TOOLS_H
#define CGAL_INTERNAL_STATIC_FILTERS_TOOLS_H

#include <CGAL/basic.h>

namespace CGAL {

template < typename ET >
class Lazy_exact_nt;

namespace internal {

// Utility function to check a posteriori that a subtraction was performed
// without rounding error.
inline bool diff_was_exact(double a, double b, double ab)
{
    return ab+b == a && a-ab == b;
}

// Auxiliary function to check if static filters can be applied, that is,
// if to_double() does not add roundoff errors.
// TODO :
// - generalize it to other number types.
// - promote it as a number type requirement ?
// - naming : is_representable_in_double() ?
//            is_representable<T>() for representable in T ?

// Current semantics :  bool fit_in_double(const NT& n, double &)
//
// - returns true means that "n" is exactly representable by a double,
//   _and_ then "returns" it in the reference.
// - it is fine to return false conservatively.

template < typename T >
inline bool fit_in_double(const T&, double&) { return false; }

inline bool fit_in_double(const double& d, double& r) { r = d; return true; }

inline bool fit_in_double(const float& f, double& r) { r = f; return true; }

inline bool fit_in_double(const int& i, double& r) { r = i; return true; }

template < typename ET >
inline bool fit_in_double(const Lazy_exact_nt<ET>&, double&);

} } // namespace CGAL::internal

#endif // CGAL_INTERNAL_STATIC_FILTERS_TOOLS_H
