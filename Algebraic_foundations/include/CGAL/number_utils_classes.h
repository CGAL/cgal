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
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

// to be included by number_utils.h

#ifndef CGAL_NUMBER_UTILS_CLASSES_H
#define CGAL_NUMBER_UTILS_CLASSES_H 1

#include <CGAL/number_type_basic.h>
#include <algorithm>
#include <utility>

CGAL_BEGIN_NAMESPACE

/* Defines functors:
- Is_zero
- Is_one
- Is_negative
- Is_positive
- Sgn
- Abs
- Compare
- Square
- Sqrt
- Div
- Gcd
- To_double
- To_interval
*/
template < class NT > struct Is_zero;
template < class NT > struct Is_one;
template < class NT > struct Is_negative;
template < class NT > struct Is_positive;
template < class NT > struct Sgn;
template < class NT > struct Abs;
template < class NT > struct Compare;
template < class NT > struct Square;
template < class NT > struct Sqrt;
template < class NT > struct Div;
template < class NT > struct Gcd;
template < class NT > struct To_double;
template < class NT > struct To_interval;


// TODO: All the following functors have to be removed since we allready have
//        them in the Real_embeddable_traits. 
template < class NT >
struct Is_zero :public Unary_function< NT, bool > {
    bool operator()( const NT& x) const
    { return CGAL_NTS is_zero( x); }
};

template < class NT >
struct Is_one :public Unary_function< NT, bool > {
    bool operator()( const NT& x) const
    { return CGAL_NTS is_one( x); }
};

template < class NT >
struct Is_negative :public Unary_function< NT, bool > {
    bool operator()( const NT& x) const
    { return CGAL_NTS is_negative( x); }
};

template < class NT >
struct Is_positive :public Unary_function< NT, bool > {
    bool operator()( const NT& x) const
    { return CGAL_NTS is_positive( x); }
};


// Sign would result in a name clash with enum.h
template < class NT >
struct Sgn : Real_embeddable_traits<NT>::Sign {};

template < class NT >
struct Abs :public Unary_function< NT, NT > {
    NT operator()( const NT& x) const
    { return CGAL_NTS abs( x); }
};

template <class NT, class Compare> struct Compare_base: public Compare {};
template <class NT> struct Compare_base<NT,Null_functor>
    :public Binary_function< NT, NT, Comparison_result > {
    Comparison_result operator()( const NT& x, const NT& y) const
    {
        if (x < y) return SMALLER;
        if (x > y) return LARGER;
        CGAL_postcondition(x == y);
        return EQUAL;
    }
};


template < class NT >
struct Compare
    :public Compare_base<NT, typename Real_embeddable_traits<NT>::Compare>{};

template < class NT >
struct Square : public Unary_function< NT, NT > {
    NT operator()( const NT& x) const
    { return CGAL_NTS square( x ); }
};

template < class NT >
struct Sqrt : public Unary_function< NT, NT > {
    NT operator()( const NT& x) const
    { return CGAL_NTS sqrt( x ); }
};

template < class NT >
struct Div : public Binary_function< NT, NT, NT > {
    NT operator()( const NT& x, const NT& y) const
    { return CGAL_NTS div( x, y ); }
};

template < class NT >
struct Gcd : public Binary_function< NT, NT, NT > {
    NT operator()( const NT& x, const NT& y) const
    { return CGAL_NTS gcd( x, y ); }
};

template < class NT >
struct To_double : public Unary_function< NT, double > {
    double operator()( const NT& x) const
    { return CGAL_NTS to_double( x ); }
};

template < class NT >
struct To_interval
    : public Unary_function< NT, std::pair<double, double> >
{
    std::pair<double, double> operator()( const NT& x) const
    { return CGAL_NTS to_interval( x ); }
};


CGAL_END_NAMESPACE

#endif // CGAL_NUMBER_UTILS_CLASSES_H
