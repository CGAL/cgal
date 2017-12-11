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
// Author(s) : Laurent Saboret and Nader Salman and Pierre Alliez

#ifndef CGAL_REMOVE_OUTLIERS_H
#define CGAL_REMOVE_OUTLIERS_H

#include <CGAL/license/Point_set_processing_3.h>


#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

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
template < typename Kernel,
           typename Tree >
typename Kernel::FT
compute_avg_knn_sq_distance_3(
    const typename Kernel::Point_3& query, ///< 3D point to project
    Tree& tree,                            ///< KD-tree
    unsigned int k)                        ///< number of neighbors
{
    // geometric types
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;

    // types for K nearest neighbors search
    typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
    typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::iterator Search_iterator;

    // Gather set of (k+1) neighboring points.
    // Perform k+1 queries (if in point set, the query point is
    // output first). Search may be aborted if k is greater
    // than number of input points.
    std::vector<Point> points; points.reserve(k+1);
    Neighbor_search search(tree,query,k+1);
    Search_iterator search_iterator = search.begin();
    unsigned int i;
    for(i=0;i<(k+1);i++)
    {
        if(search_iterator == search.end())
            break; // premature ending
        points.push_back(search_iterator->first);
        search_iterator++;
    }
    CGAL_point_set_processing_precondition(points.size() >= 1);

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

/// \ingroup PkgPointSetProcessingAlgorithms
/// Removes outliers:
/// - computes average squared distance to the K nearest neighbors,
/// - and sorts the points in increasing order of average distance.
///
/// This method modifies the order of input points so as to pack all remaining points first,
/// and returns an iterator over the first point to remove (see erase-remove idiom).
/// For this reason it should not be called on sorted containers.
///
/// \pre `k >= 2`
///
/// @tparam InputIterator iterator over input points.
/// @tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
///        It can be omitted ifthe value type of `InputIterator` is convertible to `Point_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointMap`.
///
/// @return iterator over the first point to remove.
///
/// @note There are two thresholds that can be used:
/// `threshold_percent` and `threshold_distance`. This function
/// returns the smallest number of outliers such that at least one of
/// these threshold is fullfilled. This means that if
/// `threshold_percent=100`, only `threshold_distance` is taken into
/// account; if `threshold_distance=0` only `threshold_percent` is
/// taken into account.

// This variant requires all parameters.
template <typename PointRange,
          typename NamedParameters
>
typename PointRange::iterator
remove_outliers(
  PointRange& points,
  unsigned int k, ///< number of neighbors.
  const NamedParameters& np)
{
  using boost::choose_param;
  
  // geometric types
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
  double threshold_percent = choose_param(get_param(np, internal_np::threshold_percent), 10.);
  double threshold_distance = choose_param(get_param(np, internal_np::threshold_distance), 0.);
  
  typedef typename Kernel::FT FT;
  
  // basic geometric types
  typedef typename Kernel::Point_3 Point;

  // actual type of input points
  typedef typename std::iterator_traits<typename PointRange::iterator>::value_type Enriched_point;

  // types for K nearest neighbors search structure
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2);

  CGAL_point_set_processing_precondition(threshold_percent >= 0 && threshold_percent <= 100);

  typename PointRange::iterator it;

  // Instanciate a KD-tree search.
  // Note: We have to convert each input iterator to Point_3.
  std::vector<Point> kd_tree_points; 
  for(it = points.begin(); it != points.end(); it++)
    kd_tree_points.push_back( get(point_map, *it) );
  Tree tree(kd_tree_points.begin(), kd_tree_points.end());

  // iterate over input points and add them to multimap sorted by distance to k
  std::multimap<FT,Enriched_point> sorted_points;
  for(it = points.begin(); it != points.end(); it++)
  {
    FT sq_distance = internal::compute_avg_knn_sq_distance_3<Kernel>(
      get(point_map,*it),
      tree, k);
    sorted_points.insert( std::make_pair(sq_distance, *it) );
  }

  // Replaces [points.begin(), points.end()) range by the multimap content.
  // Returns the iterator after the (100-threshold_percent) % best points.
  typename PointRange::iterator first_point_to_remove = points.begin();
  typename PointRange::iterator dst = points.begin();
  int first_index_to_remove = int(double(sorted_points.size()) * ((100.0-threshold_percent)/100.0));
  typename std::multimap<FT,Enriched_point>::iterator src;
  int index;
  for (src = sorted_points.begin(), index = 0;
       src != sorted_points.end();
       ++src, ++index)
  {
    *dst++ = src->second;
    if (index <= first_index_to_remove ||
        src->first < threshold_distance * threshold_distance)
      first_point_to_remove = dst;
  }

  return first_point_to_remove;
}

template <typename PointRange>
typename PointRange::iterator
remove_outliers(
  PointRange& points,
  unsigned int k) ///< number of neighbors.
{
  return remove_outliers (points, k, CGAL::Point_set_processing_3::parameters::all_default(points));
}
  
/// Removes outliers:
/// - computes average squared distance to the K nearest neighbors,
/// - and sorts the points in increasing order of average distance.
///
/// This method modifies the order of input points so as to pack all remaining points first,
/// and returns an iterator over the first point to remove (see erase-remove idiom).
/// For this reason it should not be called on sorted containers.
///
/// \pre `k >= 2`
///
/// @tparam InputIterator iterator over input points.
/// @tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
///        It can be omitted ifthe value type of `InputIterator` is convertible to `Point_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointMap`.
///
/// @return iterator over the first point to remove.
///
/// @note There are two thresholds that can be used:
/// `threshold_percent` and `threshold_distance`. This function
/// returns the smallest number of outliers such that at least one of
/// these threshold is fullfilled. This means that if
/// `threshold_percent=100`, only `threshold_distance` is taken into
/// account; if `threshold_distance=0` only `threshold_percent` is
/// taken into account.

// This variant requires all parameters.
template <typename InputIterator,
          typename PointMap,
          typename Kernel
>
InputIterator
remove_outliers(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointMap point_map, ///< property map: value_type of InputIterator -> Point_3
  unsigned int k, ///< number of neighbors.
  double threshold_percent, ///< maximum percentage of points to remove.
  double threshold_distance, ///< minimum distance for a point to be
                             ///< considered as outlier (distance here is the square root of the average
                             ///< squared distance to K nearest
                             ///< neighbors)
  const Kernel& /*kernel*/) ///< geometric traits.
{
  CGAL_POINT_SET_PROCESSING_DEPRECATED_V1_API("remove_outliers()");
  CGAL::Iterator_range<InputIterator> points (first, beyond);
  return remove_outliers
    (points,
     k,
     CGAL::parameters::point_map (point_map).
     threshold_percent (threshold_percent).
     threshold_distance (threshold_distance). 
     geom_traits(Kernel()));
}
  
/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename InputIterator,
          typename PointMap
>
InputIterator
remove_outliers(
  InputIterator first, ///< iterator over the first input point
  InputIterator beyond, ///< past-the-end iterator
  PointMap point_map, ///< property map: value_type of InputIterator -> Point_3
  unsigned int k, ///< number of neighbors.
  double threshold_percent, ///< percentage of points to remove
  double threshold_distance = 0.0)  ///< minimum average squared distance to K nearest neighbors
                             ///< for a point to be removed.
{
  CGAL_POINT_SET_PROCESSING_DEPRECATED_V1_API("remove_outliers()");
  CGAL::Iterator_range<InputIterator> points (first, beyond);
  return remove_outliers
    (points,
     k,
     CGAL::parameters::point_map (point_map).
     threshold_percent (threshold_percent).
     threshold_distance (threshold_distance));
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
template <typename InputIterator
>
InputIterator
remove_outliers(
  InputIterator first, ///< iterator over the first input point
  InputIterator beyond, ///< past-the-end iterator
  unsigned int k, ///< number of neighbors.
  double threshold_percent, ///< percentage of points to remove
  double threshold_distance = 0.0)  ///< minimum average squared distance to K nearest neighbors
                             ///< for a point to be removed.
{
  CGAL_POINT_SET_PROCESSING_DEPRECATED_V1_API("remove_outliers()");
  CGAL::Iterator_range<InputIterator> points (first, beyond);
  return remove_outliers
    (points,
     k,
     CGAL::parameters::threshold_percent (threshold_percent).
     threshold_distance (threshold_distance));
}
/// @endcond


} //namespace CGAL

#endif // CGAL_REMOVE_OUTLIERS_H
