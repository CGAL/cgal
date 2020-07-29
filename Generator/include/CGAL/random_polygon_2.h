// Copyright (c) 2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_RANDOM_POLYGON_2_H
#define CGAL_RANDOM_POLYGON_2_H

#include <list>
#include <set>
#include <vector>
#include <CGAL/algorithm.h>
#include <CGAL/Random_polygon_2_sweep.h>
#include <CGAL/Kernel_traits.h>

namespace CGAL {

//
// Using the provided point generator, generates a set of n points and
// produces  a simple polygon from the unique subset of points within this
// set.
//
// Each of the p possible simple polygons for the unique point set is
// generated with probability greater than 0 but the polygons are not
// generated with uniform probability.
//
template <class PointGenerator, class OutputIterator, class Traits>
OutputIterator random_polygon_2(std::size_t n,  OutputIterator result,
                                const PointGenerator& pg, const Traits& traits)
{
   typedef typename Traits::Point_2           Point_2;
   typedef std::vector<Point_2>               Vertex_list;

   Vertex_list  vertices;

   copy_n_unique(pg, n, std::back_inserter(vertices), traits);
   CGAL_assertion(!duplicate_points(vertices.begin(), vertices.end(), traits));

#ifndef CGAL_DONT_SHUFFLE_IN_RANDOM_POLYGON_2
   CGAL::cpp98::random_shuffle(vertices.begin(), vertices.end());
#endif

   make_simple_polygon(vertices.begin(), vertices.end(), traits);

   if (orientation_2(vertices.begin(), vertices.end()) == CLOCKWISE)
      std::reverse(vertices.begin(), vertices.end());

   CGAL_assertion(is_simple_2(vertices.begin(), vertices.end()));

   return std::copy(vertices.begin(), vertices.end(), result);
}

template <class PointGenerator, class OutputIterator>
inline
OutputIterator random_polygon_2( std::size_t n,  OutputIterator result,
                                 const PointGenerator& pg )
{
   typedef typename std::iterator_traits<PointGenerator>::value_type  Point_2;
   typedef typename Kernel_traits<Point_2>::Kernel   K;
   return random_polygon_2(n, result, pg, K());
}

template <class ForwardIterator, class Traits>
bool duplicate_points(ForwardIterator first, ForwardIterator beyond,
                      const Traits& )
{
   typedef typename Traits::Point_2      Point_2;
   typedef typename Traits::Less_xy_2    Less_xy_2;
   std::set<Point_2,Less_xy_2>  point_set;
   for (; first != beyond; first++)
      if (!(point_set.insert(*first)).second) return true;
   return false;
}

template <class ForwardIterator>
bool duplicate_points(ForwardIterator first, ForwardIterator beyond)
{
   typedef typename std::iterator_traits<ForwardIterator>::value_type  Point_2;
   typedef typename Kernel_traits<Point_2>::Kernel   K;
   return duplicate_points(first, beyond, K());
}


// Copies the first n points from the input iterator to the output iterator,
// removing any duplicates.  Thus fewer than n points may be inserted into
// the output iterator.
template <class InputIterator, class Size, class OutputIterator, class Traits>
OutputIterator copy_n_unique(InputIterator first, Size n,
                             OutputIterator result,
                             const Traits& )
{
   typedef typename Traits::Point_2    Point_2;
   typedef typename Traits::Less_xy_2  Less_xy_2;

   std::set<Point_2, Less_xy_2>    sorted_point_set;
   for (Size i = 0; i < n; i++)
   {
      if (sorted_point_set.insert(*first).second)
      {
          *result = *first;
          result++;
      }
      first++;
   }
   return result;
}

template <class InputIterator, class Size, class OutputIterator>
inline
OutputIterator copy_n_unique(InputIterator first, Size n,
                             OutputIterator result)
{
   typedef typename std::iterator_traits<InputIterator>::value_type  Point_2;
   typedef typename Kernel_traits<Point_2>::Kernel  K;
   return copy_n_unique(first, n, result, K());
}

} // namespace CGAL

#endif
