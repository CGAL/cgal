// Copyright (c) 2001,2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_STATIC_FILTERS_TOOLS_H
#define CGAL_STATIC_FILTERS_TOOLS_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

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
// - should it be merged in one function with to_double() ?
//   always ?  only when true ?  what if one of the values is false ?...
// - promote it as a number type requirement ?
// - naming : is_representable_in_double() ?
//            is_representable<T>() for representable in T ?

template < typename T >
inline bool fit_in_double(const T&) { return false; }

inline bool fit_in_double(const double&) { return true; }

inline bool fit_in_double(const float&) { return true; }

inline bool fit_in_double(const int&) { return true; }

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTERS_TOOLS_H
