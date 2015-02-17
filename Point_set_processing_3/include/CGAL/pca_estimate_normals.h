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
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_PCA_ESTIMATE_NORMALS_H
#define CGAL_PCA_ESTIMATE_NORMALS_H

#include <CGAL/trace.h>
#include <CGAL/Dimension.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Memory_sizer.h>

#include <iterator>
#include <list>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL
namespace internal {


/// Estimates normal direction using linear least
/// squares fitting of a plane on the K nearest neighbors.
///
/// \pre `k >= 2`
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return Computed normal. Orientation is random.
template < typename Kernel,
           typename Tree
>
typename Kernel::Vector_3
pca_estimate_normal(const typename Kernel::Point_3& query, ///< point to compute the normal at
                    Tree& tree, ///< KD-tree
                    unsigned int k) ///< number of neighbors
{
  // basic geometric types
  typedef typename Kernel::Point_3  Point;
  typedef typename Kernel::Plane_3  Plane;

  // types for K nearest neighbors search
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // Gather set of (k+1) neighboring points.
  // Perform k+1 queries (as in point set, the query point is
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
  linear_least_squares_fitting_3(points.begin(),points.end(),plane,Dimension_tag<0>());

  // output normal vector (already normalized by PCA)
  return plane.orthogonal_vector();
}

} /* namespace internal */
/// \endcond



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/// \ingroup PkgPointSetProcessing
/// Estimates normal directions of the `[first, beyond)` range of points
/// by linear least squares fitting of a plane over the k nearest neighbors.
/// The output normals are randomly oriented.
///
/// \pre `k >= 2`
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
///        It can be omitted if the value type of  `ForwardIterator` is convertible to `Point_3<Kernel>`.
/// @tparam NormalPMap is a model of `WritablePropertyMap` with value type  `Vector_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointPMap`.

// This variant requires all parameters.
template <typename ForwardIterator,
          typename PointPMap,
          typename NormalPMap,
          typename Kernel
>
void
pca_estimate_normals(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
  NormalPMap normal_pmap, ///< property map: value_type of ForwardIterator -> Vector_3.
  unsigned int k, ///< number of neighbors.
  const Kernel& /*kernel*/) ///< geometric traits.
{
  CGAL_TRACE("Calls pca_estimate_normals()\n");

  // basic geometric types
  typedef typename Kernel::Point_3 Point;

  // Input points types
  typedef typename boost::property_traits<NormalPMap>::value_type Vector;

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

  std::size_t memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("  Creates KD-tree\n");

  ForwardIterator it;

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

  memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("  Computes normals\n");

  // iterate over input points, compute and output normal
  // vectors (already normalized)
  for(it = first; it != beyond; it++)
  {
    Vector normal = internal::pca_estimate_normal<Kernel,Tree>(      
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
      get(point_pmap,it),
#else
      get(point_pmap,*it),
#endif  
      tree,
      k);

#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    put(normal_pmap, it, normal); // normal_pmap[it] = normal
#else
    put(normal_pmap, *it, normal); // normal_pmap[it] = normal
#endif 
  }

  memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("End of pca_estimate_normals()\n");
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
template <typename ForwardIterator,
          typename PointPMap,
          typename NormalPMap
>
void
pca_estimate_normals(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
  NormalPMap normal_pmap, ///< property map: value_type of ForwardIterator -> Vector_3.
  unsigned int k) ///< number of neighbors.
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  pca_estimate_normals(
    first,beyond,
    point_pmap,
    normal_pmap,
    k,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
template <typename ForwardIterator,
          typename NormalPMap
>
void
pca_estimate_normals(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  NormalPMap normal_pmap, ///< property map: value_type of ForwardIterator -> Vector_3.
  unsigned int k) ///< number of neighbors.
{
  pca_estimate_normals(
    first,beyond,
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    make_dereference_property_map(first),
#else
    make_identity_property_map(
    typename std::iterator_traits<ForwardIterator>::value_type()),
#endif
    normal_pmap,
    k);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_PCA_ESTIMATE_NORMALS_H
