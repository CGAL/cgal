// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_BARYCENTER_H
#define CGAL_BARYCENTER_H

#include <CGAL/basic.h>
#include <CGAL/Origin.h>
#include <iterator>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/Dimension_utils.h>
#include <CGAL/Dimension.h>

// Functions to compute the point given its barycentric coordinates.
// Works in 2D and 3D (and dD ?).
// Special case for the centroid.

// TODO : Note : more numerically stable variants could be implemented as well.
// TODO : Specify a traits class concept ?  Maybe not for these computations.
// TODO : Grep for "barycenter" and "centroid" in CGAL to check existing usages.
// TODO : Add barycentric_coordinates() (to the kernel, this time).

namespace CGAL {

// This one takes an iterator range over std::pair<K::Point_[23], K::FT>
template < typename InputIterator, typename K >
typename std::iterator_traits<InputIterator>::value_type::first_type
barycenter(InputIterator begin, InputIterator end, const K & )
{
  typedef typename std::iterator_traits<InputIterator>::value_type              pair;
  typedef typename pair::second_type                                            FT;
  typedef typename pair::first_type                                             Point;
  typedef typename Access::Vector<K, typename Ambient_dimension<Point, K>::type>::type  Vector;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT norm = 0;

  while (begin != end) {
    pair p = *begin++;
    v = v + p.second * (p.first - ORIGIN);
    norm += p.second;
  }

  CGAL_assertion( norm != 0 );

  return ORIGIN + v / norm;
}

// This one takes an iterator range over K::Point_[23],
// and an iterator over K::FT.
template < typename PointInputIterator, typename WeightInputIterator,
           typename K >
typename std::iterator_traits<PointInputIterator>::value_type
barycenter(PointInputIterator begin, PointInputIterator end,
           WeightInputIterator wbegin, const K & )
{
  typedef typename std::iterator_traits<PointInputIterator>::value_type         Point;
  typedef typename std::iterator_traits<WeightInputIterator>::value_type        FT;
  typedef typename Access::Vector<K, typename Ambient_dimension<Point, K>::type>::type  Vector;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT norm = 0;

  while (begin != end) {
    FT weight = *wbegin++;
    v = v + weight * (*begin++ - ORIGIN);
    norm += weight;
  }

  CGAL_assertion( norm != 0 );

  return ORIGIN + v / norm;
}

// This one takes an iterator range over std::pair<K::Point_[23], K::FT>
// And it uses Kernel_traits<> to find out its kernel.
template < typename InputIterator >
inline
typename std::iterator_traits<InputIterator>::value_type::first_type
barycenter(InputIterator begin, InputIterator end)
{
  typedef typename std::iterator_traits<InputIterator>::value_type  pair;
  typedef typename pair::first_type                                 Point;
  typedef typename Kernel_traits<Point>::Kernel                     K;

  return CGAL::barycenter(begin, end, K());
}

// This one takes an iterator range over K::Point_[23],
// and an iterator over K::FT.
// And it uses Kernel_traits<> to find out its kernel.
// To differentiate it from the others, it takes an "int" as K parameter
template < typename PointInputIterator, typename WeightInputIterator >
inline
typename std::iterator_traits<PointInputIterator>::value_type
barycenter(PointInputIterator begin, PointInputIterator end,
           WeightInputIterator wbegin, int)
{
  typedef typename std::iterator_traits<PointInputIterator>::value_type  Point;
  typedef typename Kernel_traits<Point>::Kernel                          K;

  return CGAL::barycenter(begin, end, wbegin, K());
}

} //namespace CGAL

#endif
