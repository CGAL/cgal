// Copyright (c) 2013-06  INRIA Sophia-Antipolis (France).
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
// Author(s) : Shihao Wu, Clement Jamin, Pierre Alliez 

#ifndef CGAL_BILATERAL_SMOOTH_POINT_SET_H
#define CGAL_BILATERAL_SMOOTH_POINT_SET_H

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Point_with_normal_3.h>

#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <boost/version.hpp>
#if BOOST_VERSION >= 104000
#include <boost/property_map/property_map.hpp>
#else
#include <boost/property_map.hpp>
#endif

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
const bool is_use_parallel = true;

/// \cond SKIP_IN_MANUAL

namespace CGAL {
// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace bilateral_smooth_point_set_internal{
 // Item in the Kd-tree: position (Point_3) + index
template <typename Kernel>
class Kd_tree_element : public Point_with_normal_3<Kernel>
{
public:
  unsigned int index;

  // basic geometric types
  typedef typename CGAL::Origin Origin;
  typedef CGAL::Point_with_normal_3<Kernel> Base;

  Kd_tree_element(const Origin& o = ORIGIN, unsigned int id=0)
    : Base(o), index(id)
  {}
  Kd_tree_element(const Base& p, unsigned int id=0)
    : Base(p), index(id)
  {}
  Kd_tree_element(const Kd_tree_element& other)
    : Base(other), index(other.index)
  {}
};


// Helper class for the Kd-tree
template <typename Kernel>
class Kd_tree_gt : public Kernel
{
public:
  typedef Kd_tree_element<Kernel> Point_3;
};

template <typename Kernel>
class Kd_tree_traits : public CGAL::Search_traits_3<Kd_tree_gt<Kernel> >
{
public:
  typedef typename Kernel::Point_3 PointType;
};


/// compute bilateral projection for each point
/// according to their KNN neighborhood points
/// 
/// \pre `k >= 2`, radius > 0 , sharpness_sigma > 0 && sharpness_sigma < 90
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return 

template <typename Kernel>
CGAL::Point_with_normal_3<Kernel>
compute_denoise_projection(
  const CGAL::Point_with_normal_3<Kernel>& query,      ///< 3D point to project
  const std::vector<CGAL::Point_with_normal_3<Kernel> >& neighbor_pwns,  ///< 
  typename Kernel::FT radius,          ///< accept neighborhood radius
  typename Kernel::FT sharpness_sigma  ///< control sharpness(0-90)
)
{
  CGAL_point_set_processing_precondition(radius > 0);
  CGAL_point_set_processing_precondition(sharpness_sigma > 0
                                        && sharpness_sigma < 90);

  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef CGAL::Point_with_normal_3<Kernel> Pwn;

  FT radius2 = radius * radius;

  FT weight = (FT)0.0;
  FT iradius16 = -(FT)4.0/radius2;
  FT project_dist_sum = FT(0.0);
  FT project_weight_sum = FT(0.0);
  Vector normal_sum = CGAL::NULL_VECTOR; 

  FT cos_sigma = cos(sharpness_sigma / 180.0 * 3.1415926);
  FT sharpness_bandwidth = std::pow((CGAL::max)(1e-8,1-cos_sigma), 2);

  for (unsigned int i = 0; i < neighbor_pwns.size(); ++i)
  {
    const Point& np = neighbor_pwns[i].position();
    const Vector& nn = neighbor_pwns[i].normal();

    FT dist2 = CGAL::squared_distance(query.position(), np);
    if (dist2 < radius2)
    {
      FT theta = std::exp(dist2 * iradius16);
      FT psi = std::exp(-std::pow(1-query.normal()*nn, 2)/sharpness_bandwidth);

      weight = theta * psi;

      project_dist_sum += ((query.position() - np) * nn) * weight;
      project_weight_sum += weight;
      normal_sum = normal_sum + nn * weight;
    }
  }

  Vector update_normal = normal_sum / project_weight_sum;
  update_normal = update_normal / sqrt(update_normal.squared_length());

  Point update_point = query.position() - update_normal * 
    (project_dist_sum / project_weight_sum); 

  return Pwn(update_point, update_normal);
}

/// Computes neighbors from kdtree.
///
/// \pre `k >= 2`.
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return neighbors pwn of query point.
template < typename Kernel,
           typename Tree>
std::vector<CGAL::Point_with_normal_3<Kernel> >
compute_kdtree_neighbors(
  const CGAL::Point_with_normal_3<Kernel>& query, ///< 3D point
  Tree& tree,                            ///< KD-tree
  unsigned int k                         ///< number of neighbors         
)                       
{
  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef CGAL::Point_with_normal_3<Kernel> Pwn;

  // types for K nearest neighbors search
  typedef bilateral_smooth_point_set_internal::Kd_tree_traits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // performs k + 1 queries (if unique the query point is
  // output first). search may be aborted when k is greater
  // than number of input points
  Neighbor_search search(tree,query,k+1);
  Search_iterator search_iterator = search.begin();
  ++search_iterator;
  FT max_distance = (FT)0.0;
  unsigned int i;
  std::vector<CGAL::Point_with_normal_3<Kernel> > neighbor_pwns;
  for(i = 0; i < (k+1); ++i)
  {
    if(search_iterator == search.end())
      break; // premature ending

    Pwn pwn = search_iterator->first;
    neighbor_pwns.push_back(pwn);
    ++search_iterator;
  }

  // output average max spacing
  return neighbor_pwns;
}


/// Computes max-spacing of one query point from K nearest neighbors.
///
/// \pre `k >= 2`.
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return max spacing.
template < typename Kernel,
  typename Tree >
  typename Kernel::FT
compute_max_spacing(
  const CGAL::Point_with_normal_3<Kernel>& query, ///< 3D point
  Tree& tree,                            ///< KD-tree
  unsigned int k)                        ///< number of neighbors
{
  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef CGAL::Point_with_normal_3<Kernel> Pwn;

  // types for K nearest neighbors search
  typedef bilateral_smooth_point_set_internal::Kd_tree_traits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // performs k + 1 queries (if unique the query point is
  // output first). search may be aborted when k is greater
  // than number of input points
  Neighbor_search search(tree,query,k+1);
  Search_iterator search_iterator = search.begin();
  ++search_iterator;
  FT max_distance = (FT)0.0;
  unsigned int i;
  for(i=0 ; i<(k+1) ; ++i)
  {
    if(search_iterator == search.end())
      break; // premature ending

    Pwn pwn = search_iterator->first;
    double dist2 = CGAL::squared_distance(query.position(), pwn.position());
    max_distance = (CGAL::max)(dist2, max_distance);
    ++search_iterator;
  }

  // output max spacing
  return std::sqrt(max_distance);
}

}

template <typename Kernel>
class Pwn_updater {
  typedef typename CGAL::Point_with_normal_3<Kernel> Pwn;
  typedef typename std::vector<Pwn> Pwn_set;
  typedef typename Kernel::FT FT;

  FT sharpness_sigma;
  Pwn_set* pwn_set;
  Pwn_set* update_pwn_set;
  std::vector<Pwn_set>* pwn_neighbors_set;

public:
  Pwn_updater(
    FT s, 
    Pwn_set *in,
    Pwn_set *out, 
    std::vector<Pwn_set>* neighbors) 
  : sharpness_sigma(s), 
    pwn_set(in),
    update_pwn_set(out),
    pwn_neighbors_set(neighbors)
  {}

  void operator() ( const tbb::blocked_range<size_t>& r ) const 
  { 
    for ( size_t i = r.begin(); i != r.end(); ++i ) 
    {
      (*update_pwn_set)[i] = bilateral_smooth_point_set_internal::
        compute_denoise_projection<Kernel>
        ((*pwn_set)[i], 
        (*pwn_neighbors_set)[i], 
        0.15,
        sharpness_sigma);   
    }
  }
};

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

//=============================================================================
/// \ingroup PkgPointSetProcessing
/// 
/// A function for bilateral point set denoising(smoothing) with sharp features
/// More information see: http://web.siat.ac.cn/~huihuang/EAR/EAR_page.html
/// \pre normals must be unit vectors
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` 
///         with a value_type = Point_3<Kernel>.
///         It can be omitted if ForwardIterator value_type is convertible to 
///         Point_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///      It can be omitted and deduced automatically from PointPMap value_type.
///
/// @return average point move error.

// This variant requires all parameters.
template <typename Concurrency_tag,
          typename ForwardIterator,
          typename PointPMap,
          typename NormalPMap,
          typename Kernel>
double
bilateral_smooth_point_set(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3.
  NormalPMap normal_pmap, ///< property map ForwardIterator -> Vector_3.
  const unsigned int k, ///< number of neighbors.
  const typename Kernel::FT sharpness_sigma,  ///< control sharpness(0-90)
  const Kernel& /*kernel*/) ///< geometric traits.
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename CGAL::Point_with_normal_3<Kernel> Pwn;
  typedef typename std::vector<Pwn> Pwn_set;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;

  CGAL_point_set_processing_precondition(first != beyond);
  CGAL_point_set_processing_precondition(k > 1);

  // types for K nearest neighbors search structure
  typedef bilateral_smooth_point_set_internal::
    Kd_tree_element<Kernel> Kd_tree_element;
  typedef bilateral_smooth_point_set_internal::Kd_tree_traits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;

  // copy points and normals
  Pwn_set pwn_set;
  for(ForwardIterator it = first; it != beyond; ++it)
  {
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    const Point& p = get(point_pmap, it);
    const Vector& n = get(normal_pmap, it);
#else
    const Point& p = get(point_pmap, *it);
    const Vector& n = get(normal_pmap, *it);
#endif
    pwn_set.push_back(Pwn(p, n));
  }

  unsigned int nb_points = pwn_set.size();

  CGAL::Timer task_timer;
  task_timer.start();
  std::cout << "Initialization and guess spacing: " << std::endl;

  // initiate a KD-tree search for points
  std::vector<Kd_tree_element> treeElements;
  treeElements.reserve(pwn_set.size());
  for (unsigned int i = 0 ; i < pwn_set.size(); ++i)
  {
    const Pwn& pwn = pwn_set[i];
    treeElements.push_back(Kd_tree_element(pwn, i));
  }
  Tree tree(treeElements.begin(), treeElements.end());

  // Guess spacing
  FT guess_neighbor_radius = (FT)(std::numeric_limits<double>::max)(); 
  for(unsigned int i = 0; i < nb_points; ++i)
  {
    FT max_spacing = bilateral_smooth_point_set_internal::
                     compute_max_spacing<Kernel,Tree>(pwn_set[i], tree, k);
    guess_neighbor_radius = (CGAL::max)(max_spacing, guess_neighbor_radius);
  }
  guess_neighbor_radius *= 0.95;

  CGAL::Memory_sizer::size_type memory = CGAL::Memory_sizer().virtual_size();
  std::cout << "done: " << task_timer.time() << " seconds, "
    << (memory>>20) << " Mb allocated" << std::endl;

  task_timer.start();
  std::cout << "Compute all neighbors: " << std::endl;

  // compute all neighbors
  std::vector<Pwn_set> pwn_neighbors_set;
  pwn_neighbors_set.reserve(nb_points);
  for (unsigned int i = 0 ; i < nb_points; ++i)
  {
    Pwn pwn = pwn_set[i];
    pwn_neighbors_set.push_back(bilateral_smooth_point_set_internal::
      compute_kdtree_neighbors<Kernel, Tree>(pwn, tree, k));
  }

  memory = CGAL::Memory_sizer().virtual_size();
  std::cout << "done: " << task_timer.time() << " seconds, "
    << (memory>>20) << " Mb allocated" << std::endl;
  task_timer.stop();  

  task_timer.start();
  std::cout << "Compute update points and normals: " << std::endl;
  // update points and normals
  Pwn_set update_pwn_set(nb_points);

  if(!is_use_parallel)
  {
    tbb::blocked_range<size_t> block(0, nb_points);
    Pwn_updater<Kernel> pwn_updater(sharpness_sigma,
                                    &pwn_set,
                                    &update_pwn_set,
                                    &pwn_neighbors_set);
    tbb::parallel_for(block, pwn_updater);
  }
  else
  {
    for (unsigned int i = 0 ; i < nb_points; ++i)
    {
      Pwn pwn = pwn_set[i];

      update_pwn_set[i] = bilateral_smooth_point_set_internal::
                          compute_denoise_projection<Kernel>
                          (pwn, 
                           pwn_neighbors_set[i], 
                           guess_neighbor_radius, 
                           sharpness_sigma);
    }
  }


  memory = CGAL::Memory_sizer().virtual_size();
  std::cout << "done: " << task_timer.time() << " seconds, "
    << (memory>>20) << " Mb allocated" << std::endl;
  task_timer.stop(); 

  // save results
  FT sum_move_error = 0;
  ForwardIterator it = first;
  for(unsigned int i = 0 ; it != beyond; ++it, ++i)
  {
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    Point& p = get(point_pmap, it);
    Vector& n = get(normal_pmap, it);
#else
    Point& p = get(point_pmap, *it);
    Vector& n = get(normal_pmap, *it);
#endif
    sum_move_error += CGAL::squared_distance(p, update_pwn_set[i].position());
    p = update_pwn_set[i].position();
    n = update_pwn_set[i].normal();
  }

  return sum_move_error / nb_points;
}


/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
template <typename Concurrency_tag,
          typename ForwardIterator,
          typename PointPMap,
          typename NormalPMap>
double
bilateral_smooth_point_set(
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond, ///< past-the-end input point.
  PointPMap point_pmap, ///< property map OutputIterator -> Point_3.
  NormalPMap normal_pmap, ///< property map ForwardIterator -> Vector_3.
  const unsigned int k, ///< number of neighbors.
  double sharpness_sigma  ///< control sharpness(0-90)
) ///< property map OutputIterator -> Vector_3.
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return bilateral_smooth_point_set<Concurrency_tag>(
    first, beyond,
    point_pmap,
    normal_pmap,
    k,
    sharpness_sigma,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Dereference_property_map.
template <typename Concurrency_tag,
          typename ForwardIterator,
          typename NormalPMap>
double
bilateral_smooth_point_set(
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond, ///< past-the-end input point.
  const unsigned int k, ///< number of neighbors.
  double sharpness_sigma,  ///< control sharpness(0-90)
  NormalPMap normal_pmap) ///< property map OutputIterator -> Vector_3.
{
  return bilateral_smooth_point_set<Concurrency_tag>(
    first, beyond,
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    make_dereference_property_map(first),
#else
    make_identity_property_map(
    typename std::iterator_traits<ForwardIterator>::value_type()),
#endif
    normal_pmap, 
    k,
    sharpness_sigma);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_BILATERAL_SMOOTH_POINT_SET_H
