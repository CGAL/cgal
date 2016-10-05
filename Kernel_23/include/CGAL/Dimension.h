// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// Aviv University (Israel).  All rights reserved.
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_DIMENSION_H
#define CGAL_DIMENSION_H

#include <CGAL/config.h>
#include <CGAL/Kernel_traits.h>
#include <climits>
#include <limits>
#ifdef CGAL_EIGEN3_ENABLED
#include <Eigen/Core>
#endif

namespace CGAL {

#ifdef CGAL_EIGEN3_ENABLED
const int UNKNOWN_DIMENSION=Eigen::Dynamic;
#elif defined CGAL_CXX11
const int UNKNOWN_DIMENSION=std::numeric_limits<int>::max();
#else
const int UNKNOWN_DIMENSION=(unsigned)-1/2;
#endif

// Check that dimension d1 is fine for a kernel of dimension d2.
// If d2 is unknown, any d1 is fine.
inline bool check_dimension_eq(int d1, int d2){
	return d2==UNKNOWN_DIMENSION || d1==d2;
}

// These tag classes help dispatching functions based on a geometric dimension.

template < int dim >
struct Dimension_tag
{
  static const int value = dim;
};

struct Dynamic_dimension_tag {};


namespace internal {

  template < typename D >
  struct Dim_value {
    static const int value = D::value;
  };

  template <>
  struct Dim_value <Dynamic_dimension_tag> {};

} // namespace internal


// Ambient_dimension gives access to the dimension of the ambient space of an object.

template < typename T, typename K = typename Kernel_traits<T>::Kernel >
struct Ambient_dimension
  : public internal::Dim_value< typename K::template Ambient_dimension<T>::type >
{
  typedef typename K::template Ambient_dimension<T>::type type;
};


// Feature_dimension gives access to the dimension of an object.

template < typename T, typename K = typename Kernel_traits<T>::Kernel >
struct Feature_dimension
  : public internal::Dim_value< typename K::template Feature_dimension<T>::type >
{
  typedef typename K::template Feature_dimension<T>::type type;
};

// Change the dimension
template<class D,int i=1>struct Increment_dimension {
        typedef Dynamic_dimension_tag type;
};
template<int d,int i>struct Increment_dimension<Dimension_tag<d>,i> {
        typedef Dimension_tag<d+i> type;
};

template<class D1,class D2>struct Product_dimension {
        typedef Dynamic_dimension_tag type;
};
template<int d1,int d2>struct Product_dimension<Dimension_tag<d1>,Dimension_tag<d2> > {
        typedef Dimension_tag<d1*d2> type;
};

#ifdef CGAL_EIGEN3_ENABLED
// Convert to Eigen's notion of dimension
template <class Dim_> struct Eigen_dimension {
	enum { value=Eigen::Dynamic };
};
template <int d> struct Eigen_dimension<Dimension_tag<d> > {
	enum { value=d };
};

// and convert back
template <int d> struct Dimension_eigen {
	typedef Dimension_tag<d> type;
};
template <> struct Dimension_eigen<Eigen::Dynamic> {
	typedef Dynamic_dimension_tag type;
};
#endif

} //namespace CGAL

#endif // CGAL_DIMENSION_H
