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

#ifndef CGAL_IMPROVED_JET_SMOOTH_POINT_SET_H
#define CGAL_IMPROVED_JET_SMOOTH_POINT_SET_H

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

#include <iterator>
#include <list>

namespace CGAL {

/// \cond SKIP_IN_MANUAL

// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace improved_jet_smoothing_i {


// Item in the Kd-tree: position (Point_3) + index
template <typename Kernel>
class KdTreeElement : public Kernel::Point_3
{
public:
  unsigned int index;

  // basic geometric types
  typedef typename CGAL::Origin Origin;
  typedef typename Kernel::Point_3 Point;

  KdTreeElement(const Origin& o = ORIGIN, unsigned int id=0)
    : Point(o), index(id)
  {}
  KdTreeElement(const Point& p, unsigned int id=0)
    : Point(p), index(id)
  {}
  KdTreeElement(const KdTreeElement& other)
    : Point(other), index(other.index)
  {}
};

// Helper class for the Kd-tree
template <typename Kernel>
class KdTreeGT : public Kernel
{
public:
  typedef KdTreeElement<Kernel> Point_3;
};

template <typename Kernel>
class KdTreeTraits : public CGAL::Search_traits_3<KdTreeGT<Kernel> >
{
public:
    typedef typename Kernel::Point_3 PointType;
};

/* Usage:
  typedef improved_jet_smoothing_i::KdTreeElement<Kernel> KdTreeElement;
  typedef improved_jet_smoothing_i::KdTreeTraits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;
*/


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
  const unsigned int k, // nb neighbors
  const unsigned int degree_fitting = 2,
  const unsigned int degree_monge = 2)
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;

  // types for K nearest neighbors search
  typedef improved_jet_smoothing_i::KdTreeTraits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
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
improved_jet_smooth_point(
  const typename Kernel::Point_3& pi, ///< 3D point to smooth
  const typename Kernel::Vector_3& bi, ///< bi movement
  Tree& tree, ///< KD-tree
  const std::vector<typename Kernel::Vector_3>& b,
  const unsigned int k, // nb neighbors
  typename Kernel::FT beta)
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;

  // types for K nearest neighbors search
  //typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef improved_jet_smoothing_i::KdTreeTraits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // Gather set of k neighboring points and compute the sum of their b[] values.
  // Note: we perform k+1 queries and skip the query point which is output first.
  Vector bj_sum;
  Neighbor_search search(tree,pi,k+1);
  Search_iterator search_iterator;
  for(search_iterator = search.begin(), search_iterator++; // skip pi point
      search_iterator != search.end();
      search_iterator++ )
  {
    bj_sum = bj_sum + b[search_iterator->first.index];
  }

  return pi - (beta * bi + ((1-beta)/k)*bj_sum);
}


} /* namespace improved_jet_smoothing_i */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/// Improved Laplacian smoothing (Vollmer et al)
/// on the k nearest neighbors.
/// As this method relocates the points, it
/// should not be called on containers sorted w.r.t. point locations.
///
/// \pre `k >= 2`
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadWritePropertyMap` with value type `Point_3<Kernel>`.
///        It can be omitted if the value type of `ForwardIterator` is convertible to `Point_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointPMap`.

// This variant requires all parameters.
template <typename ForwardIterator,
          typename PointPMap,
          typename Kernel
>
void
improved_jet_smooth_point_set(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
  unsigned int k, ///< number of neighbors.
  const unsigned int iter_number, ///< number of iterations.
  const Kernel& kernel, ///< geometric traits.
  typename Kernel::FT alpha,
  typename Kernel::FT beta)
{
  // Basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;

  // types for K nearest neighbors search structure
  typedef improved_jet_smoothing_i::KdTreeElement<Kernel> KdTreeElement;
  typedef improved_jet_smoothing_i::KdTreeTraits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(first != beyond);

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2);

  unsigned int i; // point index
  ForwardIterator it; // point iterator

  // Number of input points
  int nb_points = std::distance(first, beyond);

  // Instanciate a KD-tree search
  std::vector<KdTreeElement> treeElements;
  for (it = first, i=0 ; it != beyond ; ++it,++i)
  {
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    const Point& p0 = get(point_pmap,it);
#else
    const Point& p0 = get(point_pmap,*it);
#endif
    treeElements.push_back(KdTreeElement(p0,i));
  }
  Tree tree(treeElements.begin(), treeElements.end());

  std::vector<Point>  p(nb_points); // positions at step iter_n
  std::vector<Vector> b(nb_points); // ...
  for(it = first, i=0; it != beyond; it++, ++i) {
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    p[i] = get(point_pmap,it);
#else
    p[i] = get(point_pmap,*it);
#endif
  }
      

  // loop until convergence
  for(int iter_n = 0; iter_n < iter_number ; ++iter_n)
  {
      int kk = k/iter_number+1; // adaptative k

      // Iterates over input points, computes (original) Laplacian smoothing and b[].
      for(it = first, i=0; it != beyond; it++, ++i)
      {
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
          const Point& p0  = get(point_pmap,it);
#else
          const Point& p0  = get(point_pmap,*it);
#endif
          Point np   = improved_jet_smoothing_i::jet_smooth_point<Kernel>(p0,tree,kk);
          b[i]       = alpha*(np - p0) + (1-alpha)*(np - p[i]);
          p[i]       = np;
      }

      // Iterates over input points, compute and store smoothed points.
      for(it = first, i=0; it != beyond; it++, ++i)
      {
          p[i] = improved_jet_smoothing_i::improved_jet_smooth_point<Kernel>(p[i],b[i],tree,b,kk,beta);
      }
  }

  // Iterates over input points and mutates them.
  // Implementation note: the cast to Point& allows to modify only the point's position.
  for(it = first, i=0; it != beyond; it++, ++i)
  {
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    put(point_pmap, it, p[i]);
#else
    put(point_pmap, *it, p[i]);
#endif
  }
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
template <typename ForwardIterator,
          typename PointPMap,
          typename FT
>
void
improved_jet_smooth_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3
  unsigned int k, ///< number of neighbors.
  const unsigned int iter_number,
  FT alpha,
  FT beta)
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return improved_jet_smooth_point_set(
    first,beyond,
    point_pmap,
    k,
    iter_number,
    Kernel(),
    alpha, beta);
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
template <typename ForwardIterator,
          typename FT
>
void
improved_jet_smooth_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  unsigned int k, ///< number of neighbors.
  const unsigned int iter_number,
  FT alpha,
  FT beta)
{
  return improved_jet_smooth_point_set(
    first,beyond,
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    make_dereference_property_map(first),
#else
    make_identity_property_map(
    typename std::iterator_traits<ForwardIterator>::value_type()),
#endif
    k,
    iter_number,
    alpha, beta);
}
/// @endcond

/// @endcond


} //namespace CGAL

#endif // CGAL_IMPROVED_JET_SMOOTH_POINT_SET_H
