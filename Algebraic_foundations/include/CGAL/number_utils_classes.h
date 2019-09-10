// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//               : Michael Hemmer <hemmer@mpi-inf.mpg.de>

// to be included by number_utils.h

#ifndef CGAL_NUMBER_UTILS_CLASSES_H
#define CGAL_NUMBER_UTILS_CLASSES_H 1

#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <algorithm>
#include <utility>

namespace CGAL {

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

template < class NT >
struct Is_negative : Real_embeddable_traits<NT>::Is_negative {};
template < class NT >
struct Is_positive : Real_embeddable_traits<NT>::Is_positive {};
template < class NT >
struct Abs : Real_embeddable_traits<NT>::Abs{};
template < class NT >
struct To_double : Real_embeddable_traits<NT>::To_double{};
template < class NT >
struct To_interval : Real_embeddable_traits<NT>::To_interval{};

// Sign would result in a name clash with enum.h
template < class NT >
struct Sgn : Real_embeddable_traits<NT>::Sgn {};


template < class NT >
struct Square : Algebraic_structure_traits<NT>::Square{};
template < class NT >
struct Sqrt : Algebraic_structure_traits<NT>::Sqrt {};
template < class NT >
struct Div : Algebraic_structure_traits<NT>::Div{};
template < class NT >
struct Gcd : Algebraic_structure_traits<NT>::Gcd{};
template < class NT >
struct Is_one : Algebraic_structure_traits<NT>::Is_one {};

// This is due to the fact that Is_zero may be provided by 
// Algebraic_structure_traits as well as Real_embeddable_traits
// Of course it is not possible to derive from both since this 
// would cause an ambiguity. 
namespace internal{
template <class AST_Is_zero, class RET_Is_zero>
struct Is_zero_base : AST_Is_zero {} ;
template <class RET_Is_zero>
struct Is_zero_base <CGAL::Null_functor, RET_Is_zero >: RET_Is_zero {} ;
} // namespace internal
template < class NT >
struct Is_zero : 
  internal::Is_zero_base
  <typename Algebraic_structure_traits<NT>::Is_zero,
   typename Real_embeddable_traits<NT>::Is_zero>{}; 


// This is due to the fact that CGAL::Compare is used for other 
// non-realembeddable types as well.
// In this case we try to provide a default implementation
namespace internal {
template <class NT, class Compare> struct Compare_base: public Compare {};
template <class NT> struct Compare_base<NT,Null_functor>
  :public CGAL::binary_function< NT, NT, Comparison_result > {
  Comparison_result operator()( const NT& x, const NT& y) const
  {
    if (x < y) return SMALLER;
    if (x > y) return LARGER;
    CGAL_postcondition(x == y);
    return EQUAL;   
  }
};
} // namespace internal

template < class NT >
struct Compare
  :public internal::Compare_base
  <NT,typename Real_embeddable_traits<NT>::Compare>{};



} //namespace CGAL

#endif // CGAL_NUMBER_UTILS_CLASSES_H
