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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Nader Salman and Laurent Saboret 

#ifndef CGAL_GRID_SIMPLIFY_POINT_SET_H
#define CGAL_GRID_SIMPLIFY_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>


#include <CGAL/property_map.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/unordered.h>
#include <CGAL/Iterator_range.h>
#include <boost/functional/hash.hpp>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iterator>
#include <deque>
#include <algorithm>
#include <cmath>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL

namespace internal {


// Round number to multiples of epsilon
inline double round_epsilon(double value, double epsilon)
{
  return std::floor(value / epsilon);
}
  
/// Utility class for grid_simplify_point_set(): Hash_epsilon_points_3
/// defines a 3D point hash / 2 points are equal iff they belong to
/// the same cell of a grid of cell size = epsilon.
template <class Point_3, class PointMap>
struct Hash_epsilon_points_3
{
private:

    double m_epsilon;
    PointMap point_map;
    typedef typename boost::property_traits<PointMap>::value_type Point;
public:

    Hash_epsilon_points_3 (double epsilon, PointMap p_map) 
        : m_epsilon (epsilon), point_map(p_map)
    {
        CGAL_point_set_processing_precondition(epsilon > 0);
    }

  std::size_t operator() (const Point_3& a) const
  {
    const Point& pa = get(point_map,a);
    std::size_t result = boost::hash_value(round_epsilon(pa.x(), m_epsilon));
    boost::hash_combine(result, boost::hash_value(round_epsilon(pa.y(), m_epsilon)));
    boost::hash_combine(result, boost::hash_value(round_epsilon(pa.z(), m_epsilon)));
    return result;
  }

};

/// Utility class for grid_simplify_point_set(): Hash_epsilon_points_3
/// defines a 3D point equality / 2 points are equal iff they belong
/// to the same cell of a grid of cell size = epsilon.
template <class Point_3, class PointMap>
struct Equal_epsilon_points_3
{
private:

    const double m_epsilon;
    PointMap point_map;
    typedef typename boost::property_traits<PointMap>::value_type Point;
public:

    Equal_epsilon_points_3 (const double& epsilon, PointMap p_map) 
        : m_epsilon (epsilon), point_map(p_map)
    {
        CGAL_point_set_processing_precondition(epsilon > 0);
    }

    bool operator() (const Point_3& a, const Point_3& b) const
    {
      const Point& pa = get(point_map,a);
      const Point& pb = get(point_map,b);

      double ra = round_epsilon(pa.x(), m_epsilon);
      double rb = round_epsilon(pb.x(), m_epsilon);
      if (ra != rb)
        return false;
      ra = round_epsilon(pa.y(), m_epsilon);
      rb = round_epsilon(pb.y(), m_epsilon);
      if (ra != rb)
        return false;
      ra = round_epsilon(pa.z(), m_epsilon);
      rb = round_epsilon(pb.z(), m_epsilon);
      return ra == rb;
    }
};


} /* namespace internal */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/// Utility class for grid_simplify_point_set():
/// 3D points set which allows at most 1 point per cell
/// of a grid of cell size = epsilon.
///
/// Warning:
/// This class is a container sorted wrt points position
/// => you should not modify directly the order or the position of points.

template <class Point_3, class PointMap>
class Epsilon_point_set_3
  : public cpp11::unordered_set<Point_3,
                                internal::Hash_epsilon_points_3<Point_3, PointMap>,
                                internal::Equal_epsilon_points_3<Point_3, PointMap> >
{
private:

    // superclass
    typedef cpp11::unordered_set<Point_3,
                                internal::Hash_epsilon_points_3<Point_3, PointMap>,
                                internal::Equal_epsilon_points_3<Point_3, PointMap> > Base;

public:

    Epsilon_point_set_3 (double epsilon, PointMap point_map)
        : Base(10, internal::Hash_epsilon_points_3<Point_3, PointMap>(epsilon, point_map),
               internal::Equal_epsilon_points_3<Point_3, PointMap>(epsilon, point_map))
    {
        CGAL_point_set_processing_precondition(epsilon > 0);
    }

    // default copy constructor, operator =() and destructor are fine.
};

/// \endcond

/// \ingroup PkgPointSetProcessingAlgorithms
/// Merges points which belong to the same cell of a grid of cell size = `epsilon`.
///
/// This method modifies the order of input points so as to pack all remaining points first,
/// and returns an iterator over the first point to remove (see erase-remove idiom).
/// For this reason it should not be called on sorted containers.
///
/// \pre `epsilon > 0`
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
///        It can be omitted if the value type of `ForwardIterator` is convertible to `Point_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointMap`.
///
/// @return iterator over the first point to remove.

// This variant requires all parameters.
template <typename PointRange, typename NamedParameters>
typename PointRange::iterator
grid_simplify_point_set(
  PointRange& points,
  double epsilon, ///< tolerance value when merging 3D points.
  const NamedParameters& np)
{
  using boost::choose_param;
  
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::const_type PointMap;
  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());

  // actual type of input points
  typedef typename std::iterator_traits<typename PointRange::iterator>::value_type Enriched_point;

  CGAL_point_set_processing_precondition(epsilon > 0);

  // Merges points which belong to the same cell of a grid of cell size = epsilon.
  // points_to_keep[] will contain 1 point per cell; the others will be in points_to_remove[].
  Epsilon_point_set_3<Enriched_point, PointMap> points_to_keep(epsilon, point_map);
  std::deque<Enriched_point> points_to_remove;
  for (typename PointRange::iterator it = points.begin(); it != points.end(); it++)
  {
    std::pair<typename Epsilon_point_set_3<Enriched_point, PointMap>::iterator,bool> result;
    result = points_to_keep.insert(*it);
    if (!result.second) // if not inserted
      points_to_remove.push_back(*it);
  }

  // Replaces `[first, beyond)` range by the content of points_to_keep, then points_to_remove.
  typename PointRange::iterator first_point_to_remove =
    std::copy(points_to_keep.begin(), points_to_keep.end(), points.begin());
    std::copy(points_to_remove.begin(), points_to_remove.end(), first_point_to_remove);

  return first_point_to_remove;
}

template <typename PointRange>
typename PointRange::iterator
grid_simplify_point_set(PointRange& points, double epsilon)
{
  return grid_simplify_point_set (points, epsilon, CGAL::parameters::all_default());
}
  
/// Merges points which belong to the same cell of a grid of cell size = `epsilon`.
///
/// This method modifies the order of input points so as to pack all remaining points first,
/// and returns an iterator over the first point to remove (see erase-remove idiom).
/// For this reason it should not be called on sorted containers.
///
/// \pre `epsilon > 0`
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
///        It can be omitted if the value type of `ForwardIterator` is convertible to `Point_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointMap`.
///
/// @return iterator over the first point to remove.

// This variant requires all parameters.
template <typename ForwardIterator,
          typename PointMap,
          typename Kernel>
ForwardIterator grid_simplify_point_set(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3
  double epsilon, ///< tolerance value when merging 3D points.
  const Kernel& /*kernel*/) ///< geometric traits.
{
  CGAL_POINT_SET_PROCESSING_DEPRECATED_V1_API("grid_simplify_point_set()");
  CGAL::Iterator_range<ForwardIterator> points = CGAL::make_range (first, beyond);
  return grid_simplify_point_set
    (points,
     epsilon,
     CGAL::parameters::point_map (point_map).
     geom_traits (Kernel()));
}
  
/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename ForwardIterator,
          typename PointMap
>
ForwardIterator
grid_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3
  double epsilon) ///< tolerance value when merging 3D points
{
  CGAL_POINT_SET_PROCESSING_DEPRECATED_V1_API("grid_simplify_point_set()");
  CGAL::Iterator_range<ForwardIterator> points = CGAL::make_range (first, beyond);
  return grid_simplify_point_set
    (points,
     epsilon,
     CGAL::parameters::point_map (point_map));
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
template <typename ForwardIterator
>
ForwardIterator
grid_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  double epsilon) ///< tolerance value when merging 3D points
{
  CGAL_POINT_SET_PROCESSING_DEPRECATED_V1_API("grid_simplify_point_set()");
  CGAL::Iterator_range<ForwardIterator> points = CGAL::make_range (first, beyond);
  return grid_simplify_point_set
    (points,
     epsilon);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_GRID_SIMPLIFY_POINT_SET_H
