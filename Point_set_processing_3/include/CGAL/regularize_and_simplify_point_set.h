// Copyright (c) 2013-06  INRIA Sophia-Antipolis (France).
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
// Author(s) : Shihao Wu

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
/// @return iterator of the first point to downsampled points.

// This variant requires all parameters.
template <typename ForwardIterator,
          typename PointPMap,
          typename Kernel
>
ForwardIterator
regularize_and_simplify_point_set(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double retain_percentage, ///< percentage of points to retain.
  const Kernel& /*kernel*/) ///< geometric traits.
{
  CGAL_point_set_processing_precondition(retain_percentage >= 0 && retain_percentage <= 100);

  // Random shuffle
  std::random_shuffle (first, beyond);

  // Computes original(input) and sample points size 
  std::size_t nb_points_original = std::distance(first, beyond);
  std::size_t nb_points_sample = (std::size_t)(double(nb_points_original) * (retain_percentage/100.0));
  std::size_t first_index_to_sample = nb_points_original - nb_points_sample;

  // The first point iter of original and sample points
  ForwardIterator it;// point iterator
  ForwardIterator first_original_point = first;
  ForwardIterator first_sample_point = first;
  std::advance(first_sample_point, first_index_to_sample);

  //Copy sample points
  std::vector<Point> sample_points(nb_points_sample);
  unsigned int i; // sample point index
  for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
	  sample_points[i] = get(point_pmap, it);

  //Do something for sample points...



  //Copy back modified sample points to original points for output
  for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
  {
	  Point& original_p = get(point_pmap, it);
	  const Point& sample_p = sample_points[i];
	  original_p = sample_p;
  }
 
  return first_sample_point;
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
