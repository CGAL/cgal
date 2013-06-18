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

#ifndef CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H
#define CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H

#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>

namespace CGAL {

/// \ingroup PkgPointSetProcessing
/// WLOP Algorithm...
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = Point_3<Kernel>.
///        It can be omitted if ForwardIterator value_type is convertible to Point_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from PointPMap value_type.
///
/// @return iterator over the first point to sampled points.

// This variant requires all parameters.
template <typename ForwardIterator,
          typename PointPMap,
          typename Kernel
>
ForwardIterator
regularize_and_simplify_point_set(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap /*point_pmap*/, ///< property map ForwardIterator -> Point_3
  double retain_percentage, ///< percentage of points to retain.
  const Kernel& /*kernel*/) ///< geometric traits.
{
  CGAL_point_set_processing_precondition(retain_percentage >= 0 && retain_percentage <= 100);

  // Random shuffle
  std::random_shuffle (first, beyond);

  // Computes first iterator to remove
  std::size_t nb_points = std::distance(first, beyond);
  std::size_t nb_points_remain = (std::size_t)(double(nb_points) * (retain_percentage/100.0));
  std::size_t first_index_to_copy = nb_points - nb_points_remain;

  //Copy sample points
  std::vector<Point> sample_points;
  sample_points.assign(nb_points_remain, Point());
  ForwardIterator first_point_to_copy_sample_points = first;
  std::advance(first_point_to_copy_sample_points, first_index_to_copy);
  std::copy(first_point_to_copy_sample_points, beyond, sample_points.begin());

  //Do something for sample points...

  //Copy modified sample points to original points for out put
  ForwardIterator first_point_to_copy_output = first;
  std::advance(first_point_to_copy_output, first_index_to_copy);
  std::copy(sample_points.begin(), sample_points.end(), first_point_to_copy_output);

  return first_point_to_copy_output;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename ForwardIterator,
          typename PointPMap
>
ForwardIterator
regularize_and_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double removed_percentage) ///< percentage of points to remove
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return regularize_and_simplify_point_set(
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
regularize_and_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  double removed_percentage) ///< percentage of points to remove
{
  return regularize_and_simplify_point_set(
    first,beyond,
    make_dereference_property_map(first),
    removed_percentage);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H
