// Copyright (c) 2001,2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_INTERNAL_STATIC_FILTERS_TOOLS_H
#define CGAL_INTERNAL_STATIC_FILTERS_TOOLS_H

#include <CGAL/basic.h>
#include <CGAL/function_objects.h>
#include <boost/mpl/has_xxx.hpp>

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
