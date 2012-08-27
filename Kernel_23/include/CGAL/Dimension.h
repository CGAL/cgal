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

#include <CGAL/basic.h>
#include <CGAL/Kernel_traits.h>
#include <climits>

namespace CGAL {

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

} //namespace CGAL

#endif // CGAL_DIMENSION_H
