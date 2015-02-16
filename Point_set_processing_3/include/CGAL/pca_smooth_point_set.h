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
// Author(s) : Nader Salman

#ifndef CGAL_PCA_SMOOTH_POINT_SET_H
#define CGAL_PCA_SMOOTH_POINT_SET_H

#include <CGAL/Dimension.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

#include <iterator>
#include <list>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL
namespace internal {


/// Smoothes one point position using linear least
/// squares fitting of a plane (PCA) on the K
/// nearest neighbors.
///
/// \pre `k >= 2`
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return computed point
template <typename Kernel,
          typename Tree>
typename Kernel::Point_3
pca_smooth_point(
  const typename Kernel::Point_3& query, ///< 3D point to project
  Tree& tree, ///< KD-tree
  unsigned int k) ///< number of neighbors.
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Plane_3 Plane;

  // types for K nearest neighbors search
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // Gather set of (k+1) neighboring points.
  // Performs k + 1 queries (if unique the query point is
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

  // performs plane fitting by point-based PCA
  Plane plane;
  linear_least_squares_fitting_3(points.begin(),points.end(),
                                 plane,Dimension_tag<0>());

  // output projection of query point onto the plane
  return plane.projection(query);
}

} /* namespace internal */
/// \endcond


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/// \cond SKIP_IN_MANUAL
/// \ingroup PkgPointSetProcessing
/// Smoothes the `[first, beyond)` range of points using pca on the k
/// nearest neighbors and reprojection onto the plane.
/// As this method relocates the points, it
/// should not be called on containers sorted w.r.t. point locations.
///
/// \pre `k >= 2`
///
/// @tparam InputIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadWritePropertyMap` with  value type `Point_3<Kernel>`.
///        It can be omitted if the value type of `InputIterator` is convertible to `Point_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointPMap`.

// This variant requires all parameters.
template <typename InputIterator,
          typename PointPMap,
          typename Kernel
>
void
pca_smooth_point_set(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map: value_type of InputIterator -> Point_3.
  unsigned int k, ///< number of neighbors.
  const Kernel& kernel) ///< geometric traits.
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;

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

  InputIterator it;

  // Instanciate a KD-tree search.
  // Note: We have to convert each input iterator to Point_3.
  std::vector<Point> kd_tree_points; 
  for(it = first; it != beyond; it++)
  {
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    Point point = get(point_pmap, it);
#else
    Point point = get(point_pmap, *it);
#endif  
    kd_tree_points.push_back(point);
  }
  Tree tree(kd_tree_points.begin(), kd_tree_points.end());

  // Iterates over input points and mutates them.
  // Implementation note: the cast to Point& allows to modify only the point's position.
  for(it = first; it != beyond; it++)
  {
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    const Point& p = get(point_pmap, it);
    put(point_pmap, it, internal::pca_smooth_point<Kernel>(p,tree,k) );
#else
    const Point& p = get(point_pmap, *it);
    put(point_pmap, *it, internal::pca_smooth_point<Kernel>(p,tree,k) );
#endif 
  }
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
template <typename InputIterator,
          typename PointPMap
>
void
pca_smooth_point_set(
  InputIterator first, ///< iterator over the first input point
  InputIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map: value_type of InputIterator -> Point_3
  unsigned int k) ///< number of neighbors.
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  pca_smooth_point_set(
    first,beyond,
    point_pmap,
    k,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
template <typename InputIterator
>
void
pca_smooth_point_set(
  InputIterator first, ///< iterator over the first input point
  InputIterator beyond, ///< past-the-end iterator
  unsigned int k) ///< number of neighbors.
{
  pca_smooth_point_set(
    first,beyond,
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    make_dereference_property_map(first),
#else
    make_identity_property_map(
    typename std::iterator_traits<InputIterator>::value_type()),
#endif
    k);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_PCA_SMOOTH_POINT_SET_H
