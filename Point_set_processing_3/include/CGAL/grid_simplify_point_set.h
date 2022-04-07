// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Nader Salman and Laurent Saboret

#ifndef CGAL_GRID_SIMPLIFY_POINT_SET_H
#define CGAL_GRID_SIMPLIFY_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/property_map.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Iterator_range.h>
#include <functional>
#include <boost/functional/hash.hpp>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iterator>
#include <deque>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <unordered_map>

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

template <typename Point_3, typename PointMap, typename UseMap>
using Epsilon_point_set_3_base
= typename std::conditional
  <UseMap::value,
   std::unordered_map<Point_3, std::size_t,
                      internal::Hash_epsilon_points_3<Point_3, PointMap>,
                      internal::Equal_epsilon_points_3<Point_3, PointMap> >,
   std::unordered_set<Point_3,
                      internal::Hash_epsilon_points_3<Point_3, PointMap>,
                      internal::Equal_epsilon_points_3<Point_3, PointMap> > >::type;

/// Utility class for grid_simplify_point_set():
/// 3D points set which allows at most 1 (or N) point per cell
/// of a grid of cell size = epsilon.
///
/// Warning:
/// This class is a container sorted wrt points position
/// => you should not modify directly the order or the position of points.

template <class Point_3, class PointMap, class UseMap>
class Epsilon_point_set_3 : public internal::Epsilon_point_set_3_base<Point_3, PointMap, UseMap>

{
private:

  // superclass
  using Base = internal::Epsilon_point_set_3_base<Point_3, PointMap, UseMap>;

  unsigned int min_points_per_cell;

public:

  Epsilon_point_set_3 (double epsilon, PointMap point_map, unsigned int min_points_per_cell = 1)
    : Base(10, internal::Hash_epsilon_points_3<Point_3, PointMap>(epsilon, point_map),
           internal::Equal_epsilon_points_3<Point_3, PointMap>(epsilon, point_map))
    , min_points_per_cell (min_points_per_cell)
  {
    CGAL_point_set_processing_precondition(epsilon > 0);
  }

  bool insert (const Point_3& p)
  {
    return insert (p, UseMap());
  }

private:

  bool insert (const Point_3& p, const Tag_true&)
  {
    auto iter = Base::insert(std::make_pair (p, 0));
    iter.first->second ++;
    return iter.first->second == min_points_per_cell;
  }

  bool insert (const Point_3& p, const Tag_false&)
  {
    return Base::insert (p).second;
  }
};

} /* namespace internal */

/// \endcond

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessing3Algorithms
   Merges points which belong to the same cell of a grid of cell size = `epsilon`.

   This method modifies the order of input points so as to pack all remaining points first,
   and returns an iterator over the first point to remove (see erase-remove idiom).
   For this reason it should not be called on sorted containers.

   \pre `epsilon > 0`

   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range
   \param epsilon tolerance value when merging 3D points.
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadWritePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{min_points_per_cell}
       \cgalParamDescription{minimum number of points in a cell such
       that a point in this cell is kept after simplification}
       \cgalParamType{unsigned int}
       \cgalParamDefault{1}
       \cgalParamExtra{If a value greater than 1 is used, the
       algorithm also acts as an outlier filtering algorithm, by removing
       low-density areas.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \return iterator over the first point to remove.
*/
template <typename PointRange, typename NamedParameters = parameters::Default_named_parameters>
typename PointRange::iterator
grid_simplify_point_set(
  PointRange& points,
  double epsilon,
  const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef Point_set_processing_3_np_helper<PointRange, NamedParameters> NP_helper;
  typedef typename NP_helper::Point_map PointMap;
  PointMap point_map = NP_helper::get_point_map(points, np);

  unsigned int min_points_per_cell = choose_parameter(get_parameter(np, internal_np::min_points_per_cell), 1);

  // actual type of input points
  typedef typename std::iterator_traits<typename PointRange::iterator>::value_type Enriched_point;

  CGAL_point_set_processing_precondition(epsilon > 0);

  if (min_points_per_cell == 1)
  {
    // Merges points which belong to the same cell of a grid of cell size = epsilon.
    // Keep 1 point per occupied cell
    internal::Epsilon_point_set_3<Enriched_point, PointMap, Tag_false> point_set(epsilon, point_map);
    return std::partition (points.begin(), points.end(), [&](const auto& p) -> bool { return point_set.insert(p); });
  }
  // else
  // Merges points which belong to the same cell of a grid of cell size = epsilon.
  // Keep 1 point per cell occupied by at least `min_points_per_cell` points
  internal::Epsilon_point_set_3<Enriched_point, PointMap, Tag_true> point_set(epsilon, point_map, min_points_per_cell);
  return std::partition (points.begin(), points.end(), [&](const auto& p) -> bool { return point_set.insert(p); });
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_GRID_SIMPLIFY_POINT_SET_H
