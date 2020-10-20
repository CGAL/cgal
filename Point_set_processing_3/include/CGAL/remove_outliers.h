// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Laurent Saboret and Nader Salman and Pierre Alliez

#ifndef CGAL_REMOVE_OUTLIERS_H
#define CGAL_REMOVE_OUTLIERS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Point_set_processing_3/internal/Neighbor_query.h>
#include <CGAL/Point_set_processing_3/internal/Callback_wrapper.h>
#include <CGAL/for_each.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <functional>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/iterator/zip_iterator.hpp>

#include <iterator>
#include <algorithm>
#include <map>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL
namespace internal {


/// Utility function for remove_outliers():
/// Computes average squared distance to the K nearest neighbors.
///
/// \pre `k >= 2`
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return computed distance.
template <typename NeighborQuery>
typename NeighborQuery::Kernel::FT
compute_avg_knn_sq_distance_3(
  const typename NeighborQuery::Kernel::Point_3& query, ///< 3D point to project
    NeighborQuery& neighbor_query,                            ///< KD-tree
    unsigned int k,                        ///< number of neighbors
    typename NeighborQuery::Kernel::FT neighbor_radius)
{
    // geometric types
    typedef typename NeighborQuery::Kernel Kernel;
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;

    std::vector<Point> points;
    neighbor_query.get_points (query, k, neighbor_radius, std::back_inserter(points));

    // compute average squared distance
    typename Kernel::Compute_squared_distance_3 sqd;
    FT sq_distance = (FT)0.0;
    for(typename std::vector<Point>::iterator neighbor = points.begin(); neighbor != points.end(); neighbor++)
        sq_distance += sqd(*neighbor, query);
    sq_distance /= FT(points.size());
    return sq_distance;
}

} /* namespace internal */
/// \endcond


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessing3Algorithms
   Removes outliers:
   - computes average squared distance to the nearest neighbors,
   - and partitions the points either using a threshold on the of
     average distance or selecting a fixed percentage of points with
     the highest average distances

   This method modifies the order of input points so as to pack all remaining points first,
   and returns an iterator over the first point to remove (see erase-remove idiom).
   For this reason it should not be called on sorted containers.

   \pre `k >= 2`

   \tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are `Sequential_tag`,
                          `Parallel_tag`, and `Parallel_if_available_tag`.
   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param k number of neighbors
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{neighbor_radius}
       \cgalParamDescription{the spherical neighborhood radius}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{`0` (no limit)}
       \cgalParamExtra{If provided, the neighborhood of a query point is computed with a fixed spherical
                       radius instead of a fixed number of neighbors. In that case, the parameter
                       `k` is used as a limit on the number of points returned by each spherical
                       query (to avoid overly large number of points in high density areas).}
     \cgalParamNEnd

     \cgalParamNBegin{threshold_percent}
       \cgalParamDescription{the maximum percentage of points to remove}
       \cgalParamType{double}
       \cgalParamDefault{`10`}
     \cgalParamNEnd

     \cgalParamNBegin{threshold_distance}
       \cgalParamDescription{the minimum distance for a point to be considered as outlier}
       \cgalParamType{double}
       \cgalParamDefault{`0`}
       \cgalParamExtra{Distance here is the square root of the average squared distance to K-nearest neighbors}
     \cgalParamNEnd

     \cgalParamNBegin{callback}
       \cgalParamDescription{a mechanism to get feedback on the advancement of the algorithm
                             while it's running and to interrupt it if needed}
       \cgalParamType{an instance of `std::function<bool(double)>`.}
       \cgalParamDefault{unused}
       \cgalParamExtra{It is called regularly when the
                       algorithm is running: the current advancement (between 0. and 1.)
                       is passed as parameter. If it returns `true`, then the
                       algorithm continues its execution normally; if it returns
                       `false`, the algorithm is stopped, all points are left unchanged
                       and the function return `points.size()`.}
       \cgalParamExtra{The callback will be copied and therefore needs to be lightweight.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd\cgalNamedParamsEnd

   \return iterator over the first point to remove.

   \note There are two thresholds that can be used:
   `threshold_percent` and `threshold_distance`. This function
   returns the smallest number of outliers such that at least one of
   these threshold is fulfilled. This means that if
   `threshold_percent=100`, only `threshold_distance` is taken into
   account; if `threshold_distance=0` only `threshold_percent` is
   taken into account.
*/
template <typename ConcurrencyTag,
          typename PointRange,
          typename NamedParameters
>
typename PointRange::iterator
remove_outliers(
  PointRange& points,
  unsigned int k,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // geometric types
  typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  typename Kernel::FT neighbor_radius = choose_parameter(get_parameter(np, internal_np::neighbor_radius),
                                                         typename Kernel::FT(0));
  double threshold_percent = choose_parameter(get_parameter(np, internal_np::threshold_percent), 10.);
  double threshold_distance = choose_parameter(get_parameter(np, internal_np::threshold_distance), 0.);
  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                                 std::function<bool(double)>());

  typedef typename Kernel::FT FT;

  // basic geometric types
  typedef typename PointRange::iterator iterator;
  typedef typename iterator::value_type value_type;

  // types for K nearest neighbors search structure
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange&, PointMap> Neighbor_query;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2);

  CGAL_point_set_processing_precondition(threshold_percent >= 0 && threshold_percent <= 100);

  Neighbor_query neighbor_query (points, point_map);

  std::size_t nb_points = points.size();

  // iterate over input points and add them to multimap sorted by distance to k
  std::vector<std::pair<FT, value_type> > sorted_points;
  sorted_points.reserve (nb_points);
  for (const value_type& p : points)
    sorted_points.push_back(std::make_pair (FT(0), p));

  Point_set_processing_3::internal::Callback_wrapper<ConcurrencyTag>
    callback_wrapper (callback, nb_points);

  CGAL::for_each<ConcurrencyTag>
    (sorted_points,
     [&](std::pair<FT, value_type>& p) -> bool
     {
       if (callback_wrapper.interrupted())
         return false;

       p.first = internal::compute_avg_knn_sq_distance_3(
         get(point_map, p.second),
         neighbor_query, k, neighbor_radius);

       ++ callback_wrapper.advancement();
       return true;
     });

  std::size_t first_index_to_remove = std::size_t(double(sorted_points.size()) * ((100.0-threshold_percent)/100.0));

  typename std::vector<std::pair<FT, value_type> >::iterator f2r
    = sorted_points.begin();

  if (threshold_distance != FT(0))
    f2r = std::partition (sorted_points.begin(), sorted_points.end(),
                          [&threshold_distance](const std::pair<FT, value_type>& p) -> bool
                          {
                            return p.first < threshold_distance * threshold_distance;
                          });

  if (static_cast<std::size_t>(std::distance (sorted_points.begin(), f2r)) < first_index_to_remove)
  {
    std::nth_element (f2r,
                      sorted_points.begin() + first_index_to_remove,
                      sorted_points.end(),
                      [](const std::pair<FT, value_type>& v1, const std::pair<FT, value_type>& v2)
                      {
                        return v1.first<v2.first;
                      });
    f2r = sorted_points.begin() + first_index_to_remove;
  }

  // Replaces [points.begin(), points.end()) range by the sorted content.
  iterator pit = points.begin();
  iterator out = points.begin();

  for (auto sit = sorted_points.begin(); sit != sorted_points.end(); ++ sit)
  {
    *pit = sit->second;
    if (sit == f2r)
      out = pit;
    ++ pit;
  }

  callback_wrapper.join();

  // Returns the iterator on the first point to remove
  return out;
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename ConcurrencyTag, typename PointRange>
typename PointRange::iterator
remove_outliers(
  PointRange& points,
  unsigned int k) ///< number of neighbors.
{
  return remove_outliers<ConcurrencyTag> (points, k, CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_REMOVE_OUTLIERS_H
