// ======================================================================
//
// Copyright (c) 1997, 1998, 1999, 2000, 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/functional_base.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@cs.unc.edu>
//
// Basics for CGAL Functors.
// ======================================================================

#ifndef CGAL_FUNCTIONAL_BASE_H
#define CGAL_FUNCTIONAL_BASE_H

#include <functional>

CGAL_BEGIN_NAMESPACE

// +----------------------------------------------------------------------+
// | Defining a Functors Arity (== #arguments)
// +----------------------------------------------------------------------+

template < int i > struct Arity_tag { enum { arity = i }; };
// use to deduce arity of functors
// --> allows binding std functors
template < class T >
struct Arity_traits {
  typedef typename T::Arity Arity;
};

#ifndef _MSC_VER
#ifndef CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION

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

#endif // ! CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION
#endif // ! _MSC_VER

CGAL_END_NAMESPACE

#endif // CGAL_FUNCTIONAL_BASE_H
// EOF //
