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
// Author(s) : Pierre Alliez, Marc Pouget and Laurent Saboret

#ifndef CGAL_JET_SMOOTH_POINT_SET_H
#define CGAL_JET_SMOOTH_POINT_SET_H

#include <CGAL/trace.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

#include <iterator>
#include <list>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace internal {

/// \cond SKIP_IN_MANUAL

/// Smoothes one point position using jet fitting on the k
/// nearest neighbors and reprojection onto the jet.
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
jet_smooth_point(
  const typename Kernel::Point_3& query, ///< 3D point to project
  Tree& tree, ///< KD-tree
  const unsigned int k, ///< number of neighbors.
  const unsigned int degree_fitting,
  const unsigned int degree_monge)
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;

  // types for K nearest neighbors search
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // types for jet fitting
  typedef typename CGAL::Monge_via_jet_fitting<Kernel> Monge_jet_fitting;
  typedef typename Monge_jet_fitting::Monge_form Monge_form;

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

  // performs jet fitting
  Monge_jet_fitting monge_fit;
  Monge_form monge_form = monge_fit(points.begin(), points.end(),
                                    degree_fitting, degree_monge);

  // output projection of query point onto the jet
  return monge_form.origin();
}

/// \endcond

} /* namespace internal */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/// \ingroup PkgPointSetProcessing
/// Smoothes the `[first, beyond)` range of points using jet fitting on the k
/// nearest neighbors and reprojection onto the jet.
/// As this method relocates the points, it
/// should not be called on containers sorted w.r.t. point locations.
///
/// \pre `k >= 2`
///
/// @tparam InputIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = Point_3<Kernel>.
///        It can be omitted if InputIterator value_type is convertible to Point_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from PointPMap value_type.

// This variant requires all parameters.
template <typename InputIterator,
          typename PointPMap,
          typename Kernel
>
void
jet_smooth_point_set(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map InputIterator -> Point_3.
  unsigned int k, ///< number of neighbors.
  const Kernel& /*kernel*/, ///< geometric traits.
  unsigned int degree_fitting = 2, ///< fitting degree
  unsigned int degree_monge = 2)  ///< Monge degree
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
    Point point = get(point_pmap, it);
    kd_tree_points.push_back(point);
  }
  Tree tree(kd_tree_points.begin(), kd_tree_points.end());
  
  // Iterates over input points and mutates them.
  // Implementation note: the cast to Point& allows to modify only the point's position.
  for(it = first; it != beyond; it++)
  {
    Point& p = get(point_pmap, it);
    p = internal::jet_smooth_point<Kernel>(p,tree,k,degree_fitting,degree_monge);
  }
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
template <typename InputIterator,
          typename PointPMap
>
void
jet_smooth_point_set(
  InputIterator first, ///< iterator over the first input point
  InputIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map InputIterator -> Point_3
  unsigned int k, ///< number of neighbors.
  const unsigned int degree_fitting = 2,
  const unsigned int degree_monge = 2)
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  jet_smooth_point_set(
    first,beyond,
    point_pmap,
    k,
    Kernel(),
    degree_fitting,degree_monge);
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Dereference_property_map.
template <typename InputIterator
>
void
jet_smooth_point_set(
  InputIterator first, ///< iterator over the first input point
  InputIterator beyond, ///< past-the-end iterator
  unsigned int k, ///< number of neighbors.
  const unsigned int degree_fitting = 2,
  const unsigned int degree_monge = 2)
{
  jet_smooth_point_set(
    first,beyond,
    make_dereference_property_map(first),
    k,
    degree_fitting,degree_monge);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_JET_SMOOTH_POINT_SET_H

