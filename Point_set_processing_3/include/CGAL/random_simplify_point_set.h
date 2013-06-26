// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s) : Laurent Saboret

#ifndef CGAL_RANDOM_SIMPLIFY_POINT_SET_H
#define CGAL_RANDOM_SIMPLIFY_POINT_SET_H

#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>

namespace CGAL {

/// \ingroup PkgPointSetProcessing
/// Randomly deletes a user-specified fraction of the input points.
///
/// This method modifies the order of input points so as to pack all remaining points first,
/// and returns an iterator over the first point to remove (see erase-remove idiom).
/// For this reason it should not be called on sorted containers.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = Point_3<Kernel>.
///        It can be omitted if ForwardIterator value_type is convertible to Point_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from PointPMap value_type.
///
/// @return iterator over the first point to remove.

// This variant requires all parameters.
template <typename ForwardIterator,
          typename PointPMap,
          typename Kernel
>
ForwardIterator
random_simplify_point_set(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap /*point_pmap*/, ///< property map ForwardIterator -> Point_3
  double removed_percentage, ///< percentage of points to remove.
  const Kernel& /*kernel*/) ///< geometric traits.
{
  CGAL_point_set_processing_precondition(removed_percentage >= 0 && removed_percentage <= 100);

  // Random shuffle
  std::random_shuffle (first, beyond);

  // Computes first iterator to remove
  std::size_t nb_points = std::distance(first, beyond);
  std::size_t first_index_to_remove = (std::size_t)(double(nb_points) * ((100.0-removed_percentage)/100.0));
  ForwardIterator first_point_to_remove = first;
  std::advance(first_point_to_remove, first_index_to_remove);

  return first_point_to_remove;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename ForwardIterator,
          typename PointPMap
>
ForwardIterator
random_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double removed_percentage) ///< percentage of points to remove
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return random_simplify_point_set(
    first,beyond,
    point_pmap,
    removed_percentage,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Dereference_property_map.
template <typename ForwardIterator
>
ForwardIterator
random_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  double removed_percentage) ///< percentage of points to remove
{
  return random_simplify_point_set(
    first,beyond,
    make_dereference_property_map(first),
    removed_percentage);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_RANDOM_SIMPLIFY_POINT_SET_H
