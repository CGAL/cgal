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
// Author(s) : Laurent Saboret and Nader Salman and Pierre Alliez

#ifndef CGAL_REMOVE_OUTLIERS_H
#define CGAL_REMOVE_OUTLIERS_H

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

#include <iterator>
#include <algorithm>
#include <map>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace internal {

/// \cond SKIP_IN_MANUAL

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

/// \endcond

} /* namespace internal */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/// \ingroup PkgPointSetProcessing
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
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = Point_3<Kernel>.
///        It can be omitted if InputIterator value_type is convertible to Point_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from PointPMap value_type.
///
/// @return iterator over the first point to remove.

// This variant requires all parameters.
template <typename InputIterator,
          typename PointPMap,
          typename Kernel
>
InputIterator
remove_outliers(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map InputIterator -> Point_3
  unsigned int k, ///< number of neighbors.
  double threshold_percent, ///< percentage of points to remove.
  const Kernel& /*kernel*/) ///< geometric traits.
{
  // geometric types
  typedef typename Kernel::FT FT;
  
  // basic geometric types
  typedef typename Kernel::Point_3 Point;

  // actual type of input points
  typedef typename std::iterator_traits<InputIterator>::value_type Enriched_point;

  // types for K nearest neighbors search structure
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(first != beyond);

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2);

  CGAL_point_set_processing_precondition(threshold_percent >= 0 && threshold_percent <= 100);

  InputIterator it;

  // Instanciate a KD-tree search.
  // Note: We have to convert each input iterator to Point_3.
  std::vector<Point> kd_tree_points; 
  for(it = first; it != beyond; it++)
  {
    Point point = get(point_pmap, it);
    kd_tree_points.push_back(point);
  }
  Tree tree(kd_tree_points.begin(), kd_tree_points.end());

  // iterate over input points and add them to multimap sorted by distance to k
  std::multimap<FT,Enriched_point> sorted_points;
  for(it = first; it != beyond; it++)
  {
    FT sq_distance = internal::compute_avg_knn_sq_distance_3<Kernel>(get(point_pmap,it), tree, k);
    sorted_points.insert( std::make_pair(sq_distance, *it) );
  }

  // Replaces [first, beyond) range by the multimap content.
  // Returns the iterator after the (100-threshold_percent) % best points.
  InputIterator first_point_to_remove = beyond;
  InputIterator dst = first;
  int first_index_to_remove = int(double(sorted_points.size()) * ((100.0-threshold_percent)/100.0));
  typename std::multimap<FT,Enriched_point>::iterator src;
  int index;
  for (src = sorted_points.begin(), index = 0;
       src != sorted_points.end();
       ++src, ++index)
  {
    *dst++ = src->second;
    if (index == first_index_to_remove)
      first_point_to_remove = dst;
  }

  return first_point_to_remove;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename InputIterator,
          typename PointPMap
>
InputIterator
remove_outliers(
  InputIterator first, ///< iterator over the first input point
  InputIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map InputIterator -> Point_3
  unsigned int k, ///< number of neighbors.
  double threshold_percent) ///< percentage of points to remove
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return remove_outliers(
    first,beyond,
    point_pmap,
    k,threshold_percent,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Dereference_property_map.
template <typename InputIterator
>
InputIterator
remove_outliers(
  InputIterator first, ///< iterator over the first input point
  InputIterator beyond, ///< past-the-end iterator
  unsigned int k, ///< number of neighbors.
  double threshold_percent) ///< percentage of points to remove
{
  return remove_outliers(
    first,beyond,
    make_dereference_property_map(first),
    k,threshold_percent);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_REMOVE_OUTLIERS_H
