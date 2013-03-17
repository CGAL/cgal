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
// Author(s) : Nader Salman and Laurent Saboret 

#ifndef CGAL_GRID_SIMPLIFY_POINT_SET_H
#define CGAL_GRID_SIMPLIFY_POINT_SET_H

#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

#include <iterator>
#include <set>
#include <deque>
#include <algorithm>
#include <cmath>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL

namespace internal {


/// Utility class for grid_simplify_point_set():
/// Less_epsilon_points_3 defines a 3D points order / 2 points are equal
/// iff they belong to the same cell of a grid of cell size = epsilon.
template <class Point_3, class PointPMap>
struct Less_epsilon_points_3
{
private:

    double m_epsilon;
    PointPMap point_pmap;
    
public:

    Less_epsilon_points_3 (double epsilon, PointPMap p_pmap) 
        : m_epsilon (epsilon), point_pmap(p_pmap)
    {
        CGAL_point_set_processing_precondition(epsilon > 0);
    }

    bool operator() (const Point_3& a, const Point_3& b) const
    {
        typedef typename boost::property_traits<PointPMap>::value_type Point;
        
        // Round points to multiples of m_epsilon, then compare.
        Point a_n = get(point_pmap,&a),
              b_n = get(point_pmap,&b);
        
        Point rounded_a(round_epsilon(a_n.x(), m_epsilon),
                        round_epsilon(a_n.y(), m_epsilon),
                        round_epsilon(a_n.z(), m_epsilon));
        Point rounded_b(round_epsilon(b_n.x(), m_epsilon),
                        round_epsilon(b_n.y(), m_epsilon),
                        round_epsilon(b_n.z(), m_epsilon));
                        
        return (rounded_a < rounded_b);
    }

private:

    // Round number to multiples of epsilon
    static inline double round_epsilon(double value, double epsilon)
    {
        return std::floor(value/epsilon) * epsilon;
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

template <class Point_3, class PointPMap>
class Epsilon_point_set_3
  : public std::set<Point_3, internal::Less_epsilon_points_3<Point_3, PointPMap> >
{
private:

    // superclass
    typedef std::set<Point_3, internal::Less_epsilon_points_3<Point_3, PointPMap> > Base;

public:

    Epsilon_point_set_3 (double epsilon, PointPMap point_pmap)
        : Base( internal::Less_epsilon_points_3<Point_3, PointPMap>(epsilon, point_pmap) )
    {
        CGAL_point_set_processing_precondition(epsilon > 0);
    }

    // default copy constructor, operator =() and destructor are fine.
};

/// \endcond

/// \ingroup PkgPointSetProcessing
/// Merges points which belong to the same cell of a grid of cell size = `epsilon`.
///
/// This method modifies the order of input points so as to pack all remaining points first,
/// and returns an iterator over the first point to remove (see erase-remove idiom).
/// For this reason it should not be called on sorted containers.
///
/// \pre `epsilon > 0`
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
          typename Kernel>
ForwardIterator grid_simplify_point_set(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double epsilon, ///< tolerance value when merging 3D points.
  const Kernel& /*kernel*/) ///< geometric traits.
{
  // actual type of input points
  typedef typename std::iterator_traits<ForwardIterator>::value_type Enriched_point;

  CGAL_point_set_processing_precondition(epsilon > 0);

  // Merges points which belong to the same cell of a grid of cell size = epsilon.
  // points_to_keep[] will contain 1 point per cell; the others will be in points_to_remove[].
  Epsilon_point_set_3<Enriched_point, PointPMap> points_to_keep(epsilon, point_pmap);
  std::deque<Enriched_point> points_to_remove;
  for (ForwardIterator it=first ; it != beyond ; it++)
  {
      std::pair<typename Epsilon_point_set_3<Enriched_point, PointPMap>::iterator,bool> result;
      result = points_to_keep.insert(*it);
      if (!result.second) // if not inserted
        points_to_remove.push_back(*it);
  }

  // Replaces `[first, beyond)` range by the content of points_to_keep, then points_to_remove.
  ForwardIterator first_point_to_remove =
    std::copy(points_to_keep.begin(), points_to_keep.end(), first);
    std::copy(points_to_remove.begin(), points_to_remove.end(), first_point_to_remove);

  return first_point_to_remove;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename ForwardIterator,
          typename PointPMap
>
ForwardIterator
grid_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double epsilon) ///< tolerance value when merging 3D points
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return grid_simplify_point_set(
    first,beyond,
    point_pmap,
    epsilon,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Dereference_property_map.
template <typename ForwardIterator
>
ForwardIterator
grid_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  double epsilon) ///< tolerance value when merging 3D points
{
  return grid_simplify_point_set(
    first,beyond,
    make_dereference_property_map(first),
    epsilon);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_GRID_SIMPLIFY_POINT_SET_H
