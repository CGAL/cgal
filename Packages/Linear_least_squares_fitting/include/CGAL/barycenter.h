// Copyright (c) 2005  Utrecht University (The Netherlands),
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_BARYCENTER_H
#define CGAL_BARYCENTER_H

#include <CGAL/basic.h>
#include <iterator>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/Dimension_utils.h>

// Functions to compute the point given its barycentric coordinates.
// Works in 2D and 3D (and dD ?).
// Special case for the centroid.

// TODO : Note : more numerically stable variants could be implemented as well.
// TODO : Specify a traits class concept ?  Maybe not for these computations.
// TODO : Grep for "barycenter" and "centroid" in CGAL to check existing usages.
// TODO : Add barycentric_coordinates() (to the kernel, this time).

CGAL_BEGIN_NAMESPACE

// This one takes an iterator range over std::pair<K::Point_[23], K::FT>
template < typename InputIterator, typename K >
typename std::iterator_traits<InputIterator>::value_type::first_type
barycenter(InputIterator begin, InputIterator end, const K & k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type  pair;
  typedef typename pair::second_type                                FT;
  typedef typename pair::first_type                                 Point;
  typedef typename Same_dimension_vector<Point, K>::type            Vector;

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
           WeightInputIterator wbegin, const K & k)
{
  typedef typename std::iterator_traits<PointInputIterator>::value_type  Point;
  typedef typename std::iterator_traits<WeightInputIterator>::value_type FT;
  typedef typename Same_dimension_vector<Point, K>::type                 Vector;

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

CGAL_END_NAMESPACE

#endif
