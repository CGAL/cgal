// Copyright (c) 1997  Utrecht University (The Netherlands),
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
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Sylvain Pion

#ifndef CGAL_NUMBER_UTILS_FWD_H
#define CGAL_NUMBER_UTILS_FWD_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

template < class NT > struct Is_zero;
template < class NT > struct Is_one;
template < class NT > struct Is_negative;
template < class NT > struct Is_positive;
template < class NT > struct Sgn;
template < class NT > struct Abs;
template < class NT, class Compare > struct Min;
template < class NT, class Compare > struct Max;
template < class NT > struct Compare;
template < class NT > struct Square;
template < class NT > struct Sqrt;
template < class NT > struct Div;
template < class NT > struct Gcd;
template < class NT > struct To_double;
template < class NT > struct To_interval;

template <class NT>
bool is_zero(const NT& x);

template <class NT>
bool is_one(const NT& x);

template <class NT>
bool is_negative(const NT& x);

template <class NT>
bool is_positive(const NT& x);

template <class NT>
Sign sign(const NT& x);

template <class NT>
NT abs(const NT& x);

template <class NT1, class NT2>
Comparison_result compare(const NT1& n1, const NT2& n2);

template <class NT>
NT square( const NT& n);

template <class NT>
NT gcd( const NT& n1, const NT& n2);

CGAL_END_NAMESPACE

#endif // CGAL_NUMBER_UTILS_FWD_H
