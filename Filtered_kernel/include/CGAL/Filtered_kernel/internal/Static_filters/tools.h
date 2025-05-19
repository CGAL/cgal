// Copyright (c) 2001,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_INTERNAL_STATIC_FILTERS_TOOLS_H
#define CGAL_INTERNAL_STATIC_FILTERS_TOOLS_H

#include <CGAL/config.h>
#include <CGAL/function_objects.h>
#include <boost/mpl/has_xxx.hpp>

namespace CGAL {

template < typename ET >
class Lazy_exact_nt;

template <bool Protected>
class Interval_nt;

namespace internal {

// Utility function to check a posteriori that a subtraction was performed
// without rounding error.
inline bool diff_was_exact(double a, double b, double ab)
{
    return ab+b == a && a-ab == b;
}

template < typename T >
inline void init_double(double&, T* ) {}

template < typename T >
inline void init_double(double&, double&, T* ) {}

template < typename T >
inline void init_double(double&, double&, double&, T* ) {}

template < typename T >
inline void init_double(double&, double&, double&, double&, T* ) {}

template < typename T >
inline void init_double(double&, double&, double&, double&, double&, T* ) {}

template < typename T >
inline void init_double(double&, double&, double&, double&, double&, double&, T* ) {}

template < typename T >
inline void init_double(double&, double&, double&, double&, double&, double&, double&, T* ) {}

template < typename T >
inline void init_double(double&, double&, double&, double&, double&, double&, double&, double&, T* ) {}

template < typename T >
inline void init_double(double&, double&, double&,double&, double&, double&, double&, double&, double&, T* ) {}

template < typename T >
inline void init_double(double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, T* ) {}


template < typename T >
inline void init_double(double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, T* ) {}


template < typename T >
inline void init_double(double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, T* ) {}


template < typename T >
inline void init_double(double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, T* ) {}


template < typename T >
inline void init_double(double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, T* ) {}


template < typename T >
inline void init_double(double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, T* ) {}



template < typename ET >
inline void init_double(double& d0, Lazy_exact_nt<ET>* )
{d0 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, Lazy_exact_nt<ET>* )
{d0 = d1 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, double& d2, Lazy_exact_nt<ET>* )
{d0 = d1 = d2 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, double& d2, double& d3, Lazy_exact_nt<ET>* )
{d0 = d1 = d2 = d3 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, Lazy_exact_nt<ET>* )
{d0 = d1 = d2 = d3 = d4 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, Lazy_exact_nt<ET>* )
{d0 = d1 = d2 = d3 = d4 = d5 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, Lazy_exact_nt<ET>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, Lazy_exact_nt<ET>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double& d8, Lazy_exact_nt<ET>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double&d8, double& d9, Lazy_exact_nt<ET>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = d9 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double&d8, double& d9, double& d10, Lazy_exact_nt<ET>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = d9 = d10 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double&d8, double& d9, double& d10, double& d11, Lazy_exact_nt<ET>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = d9 = d10 = d11 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double&d8, double& d9, double& d10, double& d11, double& d12, Lazy_exact_nt<ET>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = d9 = d10 = d11 = d12 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double&d8, double& d9, double& d10, double& d11, double& d12, double& d13, Lazy_exact_nt<ET>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = d9 = d10 = d11 = d12 = d13 = 0;}

template < typename ET >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double&d8, double& d9, double& d10, double& d11, double& d12, double& d13, double& d14, Lazy_exact_nt<ET>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = d9 = d10 = d11 = d12 = d13 = d14 = 0;}

template < bool P >
inline void init_double(double& d0, Interval_nt<P>* )
{d0 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, Interval_nt<P>* )
{d0 = d1 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, double& d2, Interval_nt<P>* )
{d0 = d1 = d2 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, double& d2, double& d3, Interval_nt<P>* )
{d0 = d1 = d2 = d3 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, Interval_nt<P>* )
{d0 = d1 = d2 = d3 = d4 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, Interval_nt<P>* )
{d0 = d1 = d2 = d3 = d4 = d5 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, Interval_nt<P>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, Interval_nt<P>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double& d8, Interval_nt<P>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double&d8, double& d9, Interval_nt<P>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = d9 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double&d8, double& d9, double& d10, Interval_nt<P>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = d9 = d10 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double&d8, double& d9, double& d10, double& d11, Interval_nt<P>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = d9 = d10 = d11 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double&d8, double& d9, double& d10, double& d11, double& d12, Interval_nt<P>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = d9 = d10 = d11 = d12 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double&d8, double& d9, double& d10, double& d11, double& d12, double& d13, Interval_nt<P>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = d9 = d10 = d11 = d12 = d13 = 0;}

template < bool P >
inline void init_double(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5, double& d6, double& d7, double&d8, double& d9, double& d10, double& d11, double& d12, double& d13, double& d14, Interval_nt<P>* )
{d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = d8 = d9 = d10 = d11 = d12 = d13 = d14 = 0;}


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


// Auxiliary functor, to get the approximation of a kernel object:
//   - for a Point_3 of the Lazy_kernel<...>, one needs to call approx(),
//   - for a Point_3 of Simple_cartesian<double>, for example, the
//     approximation is the object itself.

namespace Static_filters_predicates {

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Kernel_object_has_approx, \
                                  Approximate_type, \
                                  false)

template <typename T, bool has_approx = Kernel_object_has_approx<T>::value>
struct Get_approx : public CGAL::Identity<T> {
  // If has_approx==false, this functor is the identity.
};

template <typename T>
struct Get_approx<T, true> {
  // If has_approx==false, this functor get .approx() on its argument.
  const typename T::Approximate_type& operator()(const T& x) const { return x.approx(); }
  typename T::Approximate_type& operator()(T& x) const { return x.approx(); }
};

} // end namespace Static_filters_predicates

} } // namespace CGAL::internal

#endif // CGAL_INTERNAL_STATIC_FILTERS_TOOLS_H
