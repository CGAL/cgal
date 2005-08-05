// Copyright (c) 1997-2004  Utrecht University (The Netherlands),
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
//                 Lutz Kettner

// Basics for CGAL Functors.

#ifndef CGAL_FUNCTIONAL_BASE_H
#define CGAL_FUNCTIONAL_BASE_H

#include <functional>

CGAL_BEGIN_NAMESPACE

// +----------------------------------------------------------------------+
// | Defining a Functors Arity (== #arguments)
// +----------------------------------------------------------------------+

template < int i > struct Arity_tag { enum { arity = i }; };

// use to deduce arity of functors --> allows binding std functors

template < class T >
struct Arity_traits {
  typedef typename T::Arity Arity;
};

// --------------------------------------------------------------------
// specializations for std functors:
//

template < class T >
struct Arity_traits< std::plus< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::minus< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::multiplies< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::divides< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::modulus< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::negate< T > > {
  typedef Arity_tag< 1 > Arity;
};
template < class T >
struct Arity_traits< std::equal_to< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::not_equal_to< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::greater< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::less< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::greater_equal< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::less_equal< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::logical_and< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::logical_or< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::logical_not< T > > {
  typedef Arity_tag< 1 > Arity;
};
template < class T >
struct Arity_traits< std::unary_negate< T > > {
  typedef Arity_tag< 1 > Arity;
};
template < class T >
struct Arity_traits< std::binary_negate< T > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T >
struct Arity_traits< std::binder1st< T > > {
  typedef Arity_tag< 1 > Arity;
};
template < class T >
struct Arity_traits< std::binder2nd< T > > {
  typedef Arity_tag< 1 > Arity;
};
template < class T1, class T2 >
struct Arity_traits< std::pointer_to_unary_function< T1, T2 > > {
  typedef Arity_tag< 1 > Arity;
};
template < class T1, class T2, class T3 >
struct Arity_traits< std::pointer_to_binary_function< T1, T2, T3 > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T1, class T2 >
struct Arity_traits< std::mem_fun_t< T1, T2 > > {
  typedef Arity_tag< 1 > Arity;
};
template < class T1, class T2, class T3 >
struct Arity_traits< std::mem_fun1_t< T1, T2, T3 > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T1, class T2 >
struct Arity_traits< std::mem_fun_ref_t< T1, T2 > > {
  typedef Arity_tag< 1 > Arity;
};
template < class T1, class T2, class T3 >
struct Arity_traits< std::mem_fun1_ref_t< T1, T2, T3 > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T1, class T2 >
struct Arity_traits< std::const_mem_fun_t< T1, T2 > > {
  typedef Arity_tag< 1 > Arity;
};
template < class T1, class T2, class T3 >
struct Arity_traits< std::const_mem_fun1_t< T1, T2, T3 > > {
  typedef Arity_tag< 2 > Arity;
};
template < class T1, class T2 >
struct Arity_traits< std::const_mem_fun_ref_t< T1, T2 > > {
  typedef Arity_tag< 1 > Arity;
};
template < class T1, class T2, class T3 >
struct Arity_traits< std::const_mem_fun1_ref_t< T1, T2, T3 > > {
  typedef Arity_tag< 2 > Arity;
};

// Versions of std::unary_function and std::binary_function which
// add the arity.
template <class Arg, class Result>
struct Unary_function : public std::unary_function<Arg, Result> {
  typedef Arity_tag< 1 > Arity;
};

template <class Arg1, class Arg2, class Result>
struct Binary_function : public std::binary_function<Arg1, Arg2, Result> {
  typedef Arity_tag< 2 > Arity;
};

CGAL_END_NAMESPACE

#endif // CGAL_FUNCTIONAL_BASE_H
