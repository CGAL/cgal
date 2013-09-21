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
// Author(s) : Shihao Wu, Clement Jamin 

#ifndef CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H
#define CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>
#include <ctime>

//for range search
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

//for AABB tree
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_point_primitive.h>
#include <CGAL/Bbox_3.h>

/// \cond SKIP_IN_MANUAL

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
  
namespace regularize_and_simplify_internal{

/// Compute average term for each sample points
/// According to their KNN neighborhood original points
/// 
/// \pre `k >= 2`, radius > 0
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return computed point

template </*typename Concurreny_tag,*/
          typename Kernel, 
          typename Tree,
          typename ForwardIterator>
typename Kernel::Vector_3
compute_average_term(
  const typename Kernel::Point_3& query, ///< 3D point to project
  Tree& aabb_tree, ///< AABB-tree
  const typename Kernel::FT radius, //accept neighborhood radius
  const std::vector<typename Kernel::FT>& density_weight_set,//if  need density
  ForwardIterator original_first_iter //
)
{
  //CGAL_point_set_processing_precondition( k > 1);
  CGAL_point_set_processing_precondition(radius > 0);
  bool is_density_weight_set_empty = density_weight_set.empty();

  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;

  // types for AABB
  typedef Kernel::Sphere_3 Circle;
  typedef std::vector<Point>::iterator Iterator;
  typedef CGAL::AABB_point_primitive<Kernel, Iterator> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> Traits_AABB;
  typedef CGAL::AABB_tree<Traits_AABB> AABB_Tree;

  std::vector<typename Primitive::Id> neighbor_original_points;
  Circle sphere_query(query, radius * radius);

  aabb_tree.all_contained_primitives(sphere_query, 
                      std::back_inserter(neighbor_original_points));

  std::cout<<"neighbor_original_points size: "<<neighbor_original_points.size()<<std::endl;

  std::vector<typename Primitive::Id>::iterator iter = neighbor_original_points.begin();
  std::vector<FT> density_set;
  std::vector<FT> dist2_set;
  FT radius2 = radius * radius;

  for (; iter != neighbor_original_points.end(); iter++)
  {
    Point& np = *(*iter);
    int original_index = std::distance(original_first_iter,
                                       *iter);

    FT dist2 = CGAL::squared_distance(query, np);
    if (!is_density_weight_set_empty)
    {
      density_set.push_back(density_weight_set[original_index]);
    }
    dist2_set.push_back(dist2);
  }

  if (neighbor_original_points.empty())
  {
    return query - CGAL::ORIGIN;
  }

  //Compute average term
  FT weight = (FT)0.0, average_weight_sum = (FT)0.0;
  FT iradius16 = -(FT)4.0 / radius2;
  Vector average = CGAL::NULL_VECTOR; 

  iter = neighbor_original_points.begin();
  int index = 0;
  for (; iter != neighbor_original_points.end(); ++iter, ++index)
  {
    Point& np = *(*iter);

    FT dist2 = dist2_set[index];
    weight = exp(dist2 * iradius16);

    if(!is_density_weight_set_empty)
    {
      weight *= density_set[index];
    }

    average_weight_sum += weight;
    average = average + (np - CGAL::ORIGIN) * weight;
  }

  // Finishing compute average term
  if (average_weight_sum > FT(1e-100))
  {
    average = average / average_weight_sum; 
  }
  else
  {
    average = query - CGAL::ORIGIN;
  }

  return average;
  // output
  //return average / average_weight_sum;
}

/// Compute repulsion term for each sample points
/// According to their KNN neighborhood sample points
/// 
/// \pre `k >= 2`, radius > 0
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return computed point
template </*typename Concurrency_tag,*/ 
          typename Kernel, 
          typename Tree,
          typename ForwardIterator>
typename Kernel::Vector_3
compute_repulsion_term(
  const typename Kernel::Point_3& query, ///< 3D point to project
  Tree& aabb_tree, ///< KD-tree
  //const unsigned int k, // nb neighbors
  const typename Kernel::FT radius, //accept neighborhood radius
  const std::vector<typename Kernel::FT>& density_weight_set, //if need density
  ForwardIterator sample_first_iter
)
{
 // CGAL_point_set_processing_precondition( k > 1);
  CGAL_point_set_processing_precondition(radius > 0);
  bool is_density_weight_set_empty = density_weight_set.empty();

  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;

  // types for AABB
  typedef Kernel::Sphere_3 Circle;
  typedef std::vector<Point>::iterator Iterator;
  typedef CGAL::AABB_point_primitive<Kernel, Iterator> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> Traits_AABB;
  typedef CGAL::AABB_tree<Traits_AABB> AABB_Tree;

  std::vector<typename Primitive::Id> neighbor_sample_points;
  Circle sphere_query(query, radius * radius);
  aabb_tree.all_contained_primitives(sphere_query, 
                            std::back_inserter(neighbor_sample_points));

  std::vector<typename Primitive::Id>::iterator iter = neighbor_sample_points.begin();
  std::vector<FT> density_set;
  std::vector<FT> dist2_set;
  FT radius2 = radius * radius;
  //parallel
  for(; iter != neighbor_sample_points.end(); iter++)
  {
    Point& np = *(*iter);
    int sample_index = std::distance(sample_first_iter, *iter);

    FT dist2 = CGAL::squared_distance(query, np);
    if (!is_density_weight_set_empty)
    {
      density_set.push_back(density_weight_set[sample_index]);
    }
    dist2_set.push_back(dist2);
  }

  if (neighbor_sample_points.empty())
  {
    return CGAL::NULL_VECTOR; 
  }

  //Compute average term
  FT weight = (FT)0.0, repulsion_weight_sum = (FT)0.0;
  FT iradius16 = -(FT)4.0/radius2;
  Vector repulsion = CGAL::NULL_VECTOR; 
  
  //parallel
  iter = neighbor_sample_points.begin();
  int index = 0;
  for(; iter != neighbor_sample_points.end(); ++iter, ++index)
  {
    Point& np = *(*iter);
    Vector diff = query - np;

    FT dist2 = dist2_set[index];
    if (dist2 < 1e-6 )  weight = 0.0;
    else
    {
      weight = std::exp(dist2 * iradius16) * 
               std::pow(FT(1.0) / std::sqrt(dist2), 2);
    }
    
    if(!is_density_weight_set_empty)
    {
      weight *= density_set[index];
    }

    repulsion_weight_sum += weight;
    repulsion = repulsion + diff * weight;
  }
  // Finishing compute repulsion term
  if (repulsion_weight_sum > FT(1e-100))
  {
    repulsion = repulsion / repulsion_weight_sum; 
  }
  else
  {
    repulsion = CGAL::NULL_VECTOR;
  }
  return repulsion;
  // output
  //return repulsion / repulsion_weight_sum;
}

/// Compute density weight for each original points,
/// according to their KNN neighborhood original points
/// 
/// \pre `k >= 2`, radius > 0
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return computed point
template </*typename Concurrency_tag, */typename Kernel, typename Tree>
typename Kernel::FT
compute_density_weight_for_original_point(
  const typename Kernel::Point_3& query, ///< 3D point to project
  Tree& aabb_tree, ///< KD-tree
  const typename Kernel::FT radius
)
{
  //CGAL_point_set_processing_precondition( k > 1);
  CGAL_point_set_processing_precondition(radius > 0);

  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;

  // types for AABB
  typedef Kernel::Sphere_3 Circle;
  typedef std::vector<Point>::iterator Iterator;
  typedef CGAL::AABB_point_primitive<Kernel, Iterator> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> Traits_AABB;
  typedef CGAL::AABB_tree<Traits_AABB> AABB_Tree;

  std::vector<typename Primitive::Id> neighbor_original_points;

  Circle sphere_query(query, radius * radius);
  aabb_tree.all_contained_primitives(sphere_query, 
                                std::back_inserter(neighbor_original_points));
  //Compute density weight
  FT radius2 = radius * radius;
  FT density_weight = (FT)1.0;
  FT iradius16 = -(FT)4.0 / radius2;

  //parallel
  std::vector<typename Primitive::Id>::iterator iter = neighbor_original_points.begin();
  for (; iter != neighbor_original_points.end(); iter++)
  {
    Point& np = *(*iter);
    FT dist2 = CGAL::squared_distance(query, np);
    density_weight += std::exp(dist2 * iradius16);
  }

  // output
  return FT(1.0) / density_weight;
}


/// Compute density weight for sample point,
/// according to their KNN neighborhood sample points
/// 
/// \pre `k >= 2`, radius > 0
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return computed point
template </*typename Concurrency_tag, */typename Kernel, typename Tree>
typename Kernel::FT
compute_density_weight_for_sample_point(
  const typename Kernel::Point_3& query, ///< 3D point to project
  Tree& aabb_tree, ///< KD-tree
  const typename Kernel::FT radius
)
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;

  // types for AABB
  typedef Kernel::Sphere_3 Circle;
  typedef std::vector<Point>::iterator Iterator;
  typedef CGAL::AABB_point_primitive<Kernel, Iterator> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> Traits_AABB;
  typedef CGAL::AABB_tree<Traits_AABB> AABB_Tree;

  std::vector<typename Primitive::Id> neighbor_sample_points;
  Circle sphere_query(query, radius * radius);
  aabb_tree.all_contained_primitives(sphere_query, 
                                   std::back_inserter(neighbor_sample_points));

  //Compute density weight
  FT radius2 = radius * radius;
  FT density_weight = (FT)1.0;
  FT iradius16 = -(FT)4.0 / radius2;
  
  //parallel
  std::vector<typename Primitive::Id>::iterator iter = neighbor_sample_points.begin();
  for (; iter != neighbor_sample_points.end(); iter++)
  {
    Point& np = *(*iter);
   
    FT dist2 = CGAL::squared_distance(query, np);
    density_weight += std::exp(dist2 * iradius16);
  }
  
  // output
  //return std::sqrt(density_weight);
  return density_weight;
}

} // namespace regularize_and_simplify_internal

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/// \ingroup PkgPointSetProcessing
/// WLOP Algorithm: The simplification algorithm can produces a set of 
/// denoised, outlier-free and evenly distributed particles over the original 
/// dense point cloud, so as to improve the reliability of other algorithms. 
///
/// The core of the algorithm is a Weighted Locally Optimal projection operator
/// with a density uniformization term. 
/// For more details, please see: http://web.siat.ac.cn/~huihuang/WLOP/WLOP_page.html
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` 
///         with a value_type = Point_3<Kernel>.
///         It can be omitted if ForwardIterator value_type is convertible to 
///         Point_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///      It can be omitted and deduced automatically from PointPMap value_type.
///
/// @return iterator of the first point to downsampled points.


// This variant requires all parameters.
template <typename ForwardIterator, typename PointPMap, typename Kernel>
ForwardIterator
regularize_and_simplify_point_set(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double retain_percentage, ///< percentage of points to retain.
  double radius, ///< number of neighbors.
  const unsigned int iter_number,///< number of iterations.
  const bool need_compute_density, ///< if needed to compute density to
                                   /// generate more rugularized result. 
  const Kernel& /*kernel*/ ///< geometric traits.
)
{
//  CGAL_point_set_processing_precondition(k > 1);
  Timer task_timer;

  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;


  // types for AABB
  typedef Kernel::Sphere_3 Circle;
  typedef std::vector<Point>::iterator Iterator;
  typedef CGAL::AABB_point_primitive<Kernel, Iterator> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> Traits_AABB;
  typedef CGAL::AABB_tree<Traits_AABB> AABB_Tree;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(first != beyond);
  CGAL_point_set_processing_precondition(retain_percentage >= 0 
                                         && retain_percentage <= 100);

  // Random shuffle
  std::random_shuffle (first, beyond);

  // Computes original(input) and sample points size 
  std::size_t nb_points_original = std::distance(first, beyond);
  std::size_t nb_points_sample = (std::size_t)(FT(nb_points_original) * 
                                 (retain_percentage/100.0));
  std::size_t first_index_to_sample = nb_points_original - nb_points_sample;

  // The first point iter of original and sample points
  ForwardIterator it;// point iterator
  ForwardIterator first_original_point = first;
  ForwardIterator first_sample_point = first;
  std::advance(first_sample_point, first_index_to_sample);

  //Copy sample points
  std::vector<Point> sample_points(nb_points_sample);
  unsigned int i; // sample point index
  //parallel
  for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
    sample_points[i] = get(point_pmap, it);

  // Initiate a KD-tree search for original points
  task_timer.start();
  
  AABB_Tree aabb_original_tree(first_original_point,
                               beyond);

  std::cout<<" build AABB_Tree Time:"<<task_timer.time()<<std::endl;
  // Compute original density weight for original points if user needed
  task_timer.reset();
  std::vector<FT> original_density_weight_set;
  if (need_compute_density)
  {
    //first resize original_density_weight_set,then parallel
    for (it = first_original_point; it != beyond ; ++it)
    {
      FT density = regularize_and_simplify_internal::
                   compute_density_weight_for_original_point<Kernel, AABB_Tree>
                                                      (get(point_pmap, it), 
                                                       aabb_original_tree, 
                                                       radius);

      original_density_weight_set.push_back(density);
    }
  }
  
  long memory = CGAL::Memory_sizer().virtual_size();
  std::cout << "compute density for original done: " << task_timer.time() << " seconds, " 
    << (memory>>20) << " Mb allocated" << std::endl;


  for (unsigned int iter_n = 0; iter_n < iter_number; iter_n++)
  {
    task_timer.reset();
    std::cout<<"AABB original tree size : "<<aabb_original_tree.size()<<std::endl;
    std::cout<<"AABB tree root node's left data"<<aabb_original_tree.bbox()<<std::endl;
    std::cout<<"Radius :"<<radius<<std::endl;
    //std::cout << "Compute average term and repulsion term " << std::endl;
    // Initiate a KD-tree search for sample points
   
    //parallel
    ForwardIterator first_sample_point = sample_points.begin();
    AABB_Tree aabb_sample_tree(sample_points.begin(),
                               sample_points.end());
    std::cout<<"build sample AABB-Tree time : "<<task_timer.time()<<std::endl;
    // Compute sample density weight for sample points if user needed
    std::vector<FT> sample_density_weight_set;
   // task_timer.start("Compute Density For Sample");
    if (need_compute_density)
    {
      //parallel
      for (i=0 ; i < sample_points.size(); i++)
      {
        FT density = regularize_and_simplify_internal::
                     compute_density_weight_for_sample_point<Kernel, AABB_Tree>
                     (sample_points[i], aabb_sample_tree, radius);

        sample_density_weight_set.push_back(density);
      }
    }
    long memory = CGAL::Memory_sizer().virtual_size();
    std::cout << "compute density for sample done: " << task_timer.time() << " seconds, " 
      << (memory>>20) << " Mb allocated" << std::endl;
    task_timer.reset();

    // Compute average term and repulsion term for each sample points ,
    // then update each sample points
    std::vector<Vector> average_set(nb_points_sample);
    std::vector<Vector> repulsion_set(nb_points_sample);
    //task_timer.start("Compute Average Term");
    //parallel
    for (i = 0; i < sample_points.size(); i++)
    {
      Point& p = sample_points[i];
      average_set[i] = regularize_and_simplify_internal::
                       compute_average_term<Kernel, AABB_Tree, ForwardIterator>
                                           (p, 
                                            aabb_original_tree, 
                                            radius, 
                                            original_density_weight_set,
                                            first_original_point);
    }

    std::cout<<"compute average term time: "<<task_timer.time()<<std::endl;
    task_timer.reset();

    //task_timer.start("Compute Repulsion Term");
    //parallel
    for (i = 0; i < sample_points.size(); i++)
    {
      Point& p = sample_points[i];
      repulsion_set[i] = regularize_and_simplify_internal::
                         compute_repulsion_term<Kernel, AABB_Tree, ForwardIterator>
                                               (p, 
                                                aabb_sample_tree, 
                                                radius, 
                                                sample_density_weight_set,
                                                first_sample_point);
    }

    for (i = 0; i < sample_points.size(); i++)
    {
      Point& p = sample_points[i];
      p = CGAL::ORIGIN + average_set[i] + (FT)0.5 * repulsion_set[i];
    }

    memory = CGAL::Memory_sizer().virtual_size();
    std::cout << "compute_repulsion_term done: " << task_timer.time() << " seconds, " 
              << (memory>>20) << " Mb allocated" << std::endl;
   
    std::cout << "iterate: " << iter_n + 1 << std::endl << std::endl;
  }

  //Copy back modified sample points to original points for output
  //parallel
  for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
  {
    Point& original_p = get(point_pmap, it);
    const Point& sample_p = sample_points[i];
    original_p = sample_p;
  }

  return first_sample_point;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename ForwardIterator, typename PointPMap>
ForwardIterator
regularize_and_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double retain_percentage, ///< percentage of points to retain.
  double radius, ///< number of neighbors.
  const unsigned int iter_number, ///< number of iterations.
  const bool need_compute_density  ///< if needed to compute density to 
                                   /// generate more rugularized result
) 
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return regularize_and_simplify_point_set(
    first,beyond,
    point_pmap,
    retain_percentage,
    radius,
    iter_number,
    need_compute_density,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Dereference_property_map
template <typename ForwardIterator>
ForwardIterator
regularize_and_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  double retain_percentage, ///< percentage of points to retain.
  double radius, ///< number of neighbors.
  const unsigned int iter_number, ///< number of iterations.
  const bool need_compute_density ///< if needed to compute density to 
                                  /// generate more rugularized result                               
) 
{
  return regularize_and_simplify_point_set(
    first,beyond,
    make_dereference_property_map(first),
    retain_percentage, radius, iter_number, need_compute_density);
}
/// @endcond

} //namespace CGAL

#endif // CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H
