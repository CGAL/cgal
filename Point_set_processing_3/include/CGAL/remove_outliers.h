// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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

CGAL_BEGIN_NAMESPACE


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace CGALi {


/// Utility function for remove_outliers():
/// Computes average squared distance to the K nearest neighbors.
///
/// @commentheading Precondition: k >= 2.
///
/// @commentheading Template Parameters:
/// @param Kernel Geometric traits class.
/// @param Tree KD-tree.
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


} /* namespace CGALi */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/// Removes outliers:
/// - computes average squared distance to the K nearest neighbors,
/// - and sorts the points in increasing order of average distance.
///
/// This method modifies the order of input points so as to pack all remaining points first,
/// and returns an iterator over the first point to remove (see erase-remove idiom).
/// For this reason it should not be called on sorted containers.
///
/// @commentheading Precondition: k >= 2.
///
/// @commentheading Template Parameters:
/// @param ForwardIterator iterator over input points.
/// @param PointPMap is a model of boost::ReadablePropertyMap with a value_type = Point_3<Kernel>.
///        It can be omitted if ForwardIterator value_type is convertible to Point_3<Kernel>.
/// @param Kernel Geometric traits class.
///        It can be omitted and deduced automatically from PointPMap value_type.
///
/// @return iterator over the first point to remove.

// This variant requires all parameters.
template <typename ForwardIterator,
          typename PointPMap,
          typename Kernel
>
ForwardIterator
remove_outliers(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  unsigned int k, ///< number of neighbors.
  double threshold_percent, ///< percentage of points to remove.
  const Kernel& kernel) ///< geometric traits.
{
  // geometric types
  typedef typename Kernel::FT FT;

  // actual type of input points
  typedef typename std::iterator_traits<ForwardIterator>::value_type Enriched_point;

  // types for K nearest neighbors search structure
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(first != beyond);

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2);

  CGAL_point_set_processing_precondition(threshold_percent >= 0 && threshold_percent <= 100);

  // instanciate a KD-tree search
  Tree tree(first,beyond);

  // iterate over input points and add them to multimap sorted by distance to k
  std::multimap<FT,Enriched_point> sorted_points;
  for(ForwardIterator it = first; it != beyond; it++)
  {
    FT sq_distance = CGALi::compute_avg_knn_sq_distance_3<Kernel>(get(point_pmap,it), tree, k);
    sorted_points.insert( std::make_pair(sq_distance, *it) );
  }

  // Replaces [first, beyond) range by the multimap content.
  // Returns the iterator after the (100-threshold_percent) % best points.
  ForwardIterator first_point_to_remove = beyond;
  ForwardIterator dst = first;
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
template <typename ForwardIterator,
          typename PointPMap
>
ForwardIterator
remove_outliers(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
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
template <typename ForwardIterator
>
ForwardIterator
remove_outliers(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  unsigned int k, ///< number of neighbors.
  double threshold_percent) ///< percentage of points to remove
{
  return remove_outliers(
    first,beyond,
    make_dereference_property_map(first),
    k,threshold_percent);
}
/// @endcond


CGAL_END_NAMESPACE

#endif // CGAL_REMOVE_OUTLIERS_H

