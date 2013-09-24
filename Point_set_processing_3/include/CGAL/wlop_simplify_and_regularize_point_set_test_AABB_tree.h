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

#ifndef CGAL_wlop_simplify_and_regularize_point_set_H
#define CGAL_wlop_simplify_and_regularize_point_set_H

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

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/tbbmalloc_proxy.h>

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
  
namespace simplify_and_regularize_internal{

/// Compute average term for each sample points
/// According to their KNN neighborhood original points
/// 
/// \pre `k >= 2`, radius > 0
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return computed point

template <typename Concurrency_tag,
          typename Kernel, 
          typename Tree,
          typename RandomAccessIterator>
typename Kernel::Vector_3
compute_average_term(
  const typename Kernel::Point_3& query, ///< 3D point to project
  Tree& aabb_tree,                       ///< AABB-tree
  const typename Kernel::FT radius,      //accept neighborhood radius
  const std::vector<typename Kernel::FT>& density_weight_set,//if  need density
  RandomAccessIterator original_first_iter
)
{
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

  std::vector<typename Primitive::Id>::iterator iter = neighbor_original_points.begin();
  std::vector<FT> density_set;
  std::vector<FT> dist2_set;
  FT radius2 = radius * radius;

  for (; iter != neighbor_original_points.end(); iter++)
  {
    Point& np = *(*iter);
    int original_index = std::distance(original_first_iter,*iter);

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
template <typename Kernel, 
          typename Tree,
          typename RandomAccessIterator>
typename Kernel::Vector_3
compute_repulsion_term(
  const typename Kernel::Point_3& query, ///< 3D point to project
  Tree& aabb_tree,                       ///<AABB-tree
  const typename Kernel::FT radius,      //accept neighborhood radius
  const std::vector<typename Kernel::FT>& density_weight_set, //if need density
  RandomAccessIterator sample_first_iter
)
{
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
template <typename Kernel, typename Tree>
typename Kernel::FT
compute_density_weight_for_original_point(
  const typename Kernel::Point_3& query, ///< 3D point to project
  Tree& aabb_tree, ///< AABB-tree
  const typename Kernel::FT radius
)
{
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
template <typename Kernel, typename Tree>
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

  std::vector<typename Primitive::Id>::iterator iter = neighbor_sample_points.begin();
  for (; iter != neighbor_sample_points.end(); iter++)
  {
    Point& np = *(*iter);
   
    FT dist2 = CGAL::squared_distance(query, np);
    density_weight += std::exp(dist2 * iradius16);
  }
  
  return density_weight;
}

} // namespace simplify_and_regularize_internal

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
/// @tparam RandomAccessIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` 
///         with a value_type = Point_3<Kernel>.
///         It can be omitted if RandomAccessIterator value_type is convertible to 
///         Point_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///      It can be omitted and deduced automatically from PointPMap value_type.
///
/// @return iterator of the first point to down sampled points.


// This variant requires all parameters.
template <typename Concurrency_tag,
          typename OutputIteratorValueType,
          typename OutputIterator,     //add output iterator
          typename RandomAccessIterator ,
          typename PointPMap,
          typename Kernel>
void
wlop_simplify_and_regularize_point_set(
  RandomAccessIterator first,  ///< iterator over the first input point.
  RandomAccessIterator beyond, ///< past-the-end iterator over the input points.
  OutputIterator output,        //add back-inserter
  PointPMap point_pmap, ///< property map RandomAccessIterator  -> Point_3
  PointPMap point_pmap_output,
  double retain_percentage, ///< percentage of points to retain.
  double radius, ///< number of neighbors.
  const unsigned int iter_number,///< number of iterations.
  const bool need_compute_density, ///< if needed to compute density to
                                   /// generate more regularized result. 
  const Kernel& /*kernel*/ ///< geometric traits.
)
{
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
  RandomAccessIterator it;// point iterator
  RandomAccessIterator first_original_point = first;
  RandomAccessIterator first_sample_point = first;
  std::advance(first_sample_point, first_index_to_sample);

  //Copy sample points
  std::vector<Point> sample_points(nb_points_sample);
  unsigned int i; // sample point index

  for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
  {
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    sample_points[i] = get(point_pmap, it);
#else
    sample_points[i] = get(point_pmap, *it);
#endif
  }

  task_timer.start();

  // Initiate a AABB_Tree search for original points
  AABB_Tree aabb_original_tree(first_original_point,
                               beyond);

  // Compute original density weight for original points if user needed
  std::vector<FT> original_density_weight_set(nb_points_original);

  if (need_compute_density)
  {
    //parallel
#ifdef CGAL_LINKED_WITH_TBB
    if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
    {
      tbb::parallel_for(
        tbb::blocked_range<size_t>(0,nb_points_original),
        [&](const tbb::blocked_range<size_t>& r)
      {
        for (size_t i = r.begin(); i< r.end(); ++i)
        {
          RandomAccessIterator cur = first;
          std::advance(cur, i);
            FT density = simplify_and_regularize_internal::
                   compute_density_weight_for_original_point<Kernel, AABB_Tree>
                                                (
                                                #ifdef CGAL_USE_PROPERTY_MAPS_API_V1
                                                  get(point_pmap, cur),
                                                #else
                                                  get(point_pmap, *cur),
                                                #endif 
                                                  aabb_original_tree, 
                                                  radius);

            original_density_weight_set[i] = density;
        }
      }
      );
    }else
#endif
    {
      for (it = first_original_point, i = 0; it != beyond ; ++it, ++i)
      {
        FT density = simplify_and_regularize_internal::
                      compute_density_weight_for_original_point<Kernel, AABB_Tree>
                                              (
                                              #ifdef CGAL_USE_PROPERTY_MAPS_API_V1
                                                get(point_pmap, it),
                                              #else
                                                get(point_pmap, *it),
                                              #endif  
                                                aabb_original_tree, 
                                                radius);

        original_density_weight_set[i] = density;
      }
    }
  }
  
  long memory = CGAL::Memory_sizer().virtual_size();
  std::cout << "compute density for original done: " << task_timer.time() << " seconds, " 
    << (memory>>20) << " Mb allocated" << std::endl << std::endl;


  for (unsigned int iter_n = 0; iter_n < iter_number; iter_n++)
  {
    task_timer.reset();
    RandomAccessIterator first_sample_point = sample_points.begin();
    AABB_Tree aabb_sample_tree(sample_points.begin(),
                               sample_points.end());

    // Compute sample density weight for sample points if user needed
    std::vector<FT> sample_density_weight_set;
    if (need_compute_density)
    {
      for (i=0 ; i < sample_points.size(); i++)
      {
        FT density = simplify_and_regularize_internal::
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

    //parallel
#ifdef CGAL_LINKED_WITH_TBB
    if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
    {
      tbb::parallel_for(
        tbb::blocked_range<size_t>(0,nb_points_sample),
        [&](const tbb::blocked_range<size_t>& r)
      {
        for (size_t i = r.begin(); i< r.end(); ++i)
        {
          Point& p = sample_points[i];
          average_set[i] = simplify_and_regularize_internal::
            compute_average_term<Concurrency_tag, Kernel, AABB_Tree, RandomAccessIterator>
            (p, 
            aabb_original_tree, 
            radius, 
            original_density_weight_set,
            first_original_point);
        }
      }
      );
    }else
#endif
    {
      for (i = 0; i < sample_points.size(); i++)
      {
        Point& p = sample_points[i];
        average_set[i] = simplify_and_regularize_internal::
          compute_average_term<Concurrency_tag, Kernel, AABB_Tree, RandomAccessIterator>
          (p, 
          aabb_original_tree, 
          radius, 
          original_density_weight_set,
          first_original_point);
      }
    }

    std::cout<<"compute average term time: "<<task_timer.time()<<std::endl;
    task_timer.reset();

    //parallel
#ifdef CGAL_LINKED_WITH_TBB
    if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
    {
      tbb::parallel_for(
        tbb::blocked_range<size_t>(0,nb_points_sample),
        [&](const tbb::blocked_range<size_t>& r)
      {
        for (size_t i = r.begin(); i< r.end(); ++i)
        {
          Point& p = sample_points[i];
          repulsion_set[i] = simplify_and_regularize_internal::
            compute_repulsion_term<Kernel, AABB_Tree, RandomAccessIterator>
            (p, 
            aabb_sample_tree, 
            radius, 
            sample_density_weight_set,
            first_sample_point);

          p = CGAL::ORIGIN + average_set[i] + (FT)0.5 * repulsion_set[i];
        }
        
      }
      );
    }else
#endif
    {
      for (i = 0; i < sample_points.size(); i++)
      {
        Point& p = sample_points[i];
        repulsion_set[i] = simplify_and_regularize_internal::
          compute_repulsion_term<Kernel, AABB_Tree, RandomAccessIterator>
          (p, 
          aabb_sample_tree, 
          radius, 
          sample_density_weight_set,
          first_sample_point);

        p = CGAL::ORIGIN + average_set[i] + (FT)0.5 * repulsion_set[i];
      }
    }

    memory = CGAL::Memory_sizer().virtual_size();
    std::cout << "compute_repulsion_term done: " << task_timer.time() << " seconds, " 
              << (memory>>20) << " Mb allocated" << std::endl;
   
    std::cout << "iterate: " << iter_n + 1 << std::endl << std::endl;
  }

  for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
  {
    const Point& sample_p = sample_points[i];
    *output++ = sample_p;
  }

  return;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename Concurrency_tag,
          typename OutputIteratorValueType,
          typename OutputIterator,     //add output iterator
          typename RandomAccessIterator, 
          typename PointPMap>
void 
wlop_simplify_and_regularize_point_set(
  RandomAccessIterator  first, ///< iterator over the first input point
  RandomAccessIterator  beyond, ///< past-the-end iterator
  OutputIterator output,        //add back-inserter
  PointPMap point_pmap,      ///< property map RandomAccessIterator  -> Point_3
  PointPMap point_pmap_output,
  double retain_percentage, ///< percentage of points to retain
  double neighbor_radius, ///< size of neighbors.
  const unsigned int max_iter_number, ///< number of iterations.
  const bool need_compute_density  ///< if needed to compute density 
                                   ///  to generate more rugularized result.                                 
) 
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  typedef typename value_type_traits<OutputIterator>::type OutputIteratorType;
  return wlop_simplify_and_regularize_point_set<Concurrency_tag, OutputIteratorType>(
    first, beyond,
    output,
    point_pmap,
    point_pmap_output,
    retain_percentage,
    neighbor_radius,
    max_iter_number,
    need_compute_density,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
/// This variant creates a default point property map=Dereference_property_map.
template <typename Concurrency_tag, 
          typename OutputIteratorValueType,
          typename OutputIterator,     //add output iterator
          typename RandomAccessIterator,
          typename PointPMap>
void
wlop_simplify_and_regularize_point_set(
  RandomAccessIterator  first,  ///< iterator over the first input point
  RandomAccessIterator  beyond, ///< past-the-end iterator
  OutputIterator output,        //add back-inserter
  PointPMap point_pmap,      ///< property map RandomAccessIterator  -> Point_3
  double retain_percentage = 5, ///< percentage of points to retain
  double neighbor_radius = 0.25, ///< size of neighbors.
  const unsigned int max_iter_number = 35, ///< number of iterations.
  const bool need_compute_density = true ///< if needed to compute density to   
                                          /// generate more uniform result. 
)
{
   return wlop_simplify_and_regularize_point_set
    <Concurrency_tag, OutputIteratorValueType>(
    first, beyond,
    output,
    point_pmap,
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    make_dereference_property_map(output),
#else
    make_identity_property_map(OutputIteratorValueType()),
#endif
    retain_percentage, 
    neighbor_radius, 
    max_iter_number, 
    need_compute_density);
}
/// @endcond


/// @cond SKIP_IN_MANUAL
/// This variant creates a default point property map=Dereference_property_map.
template <typename Concurrency_tag, 
          typename OutputIterator,     //add output iterator
          typename RandomAccessIterator >
void
wlop_simplify_and_regularize_point_set(
  RandomAccessIterator  first,  ///< iterator over the first input point
  RandomAccessIterator  beyond, ///< past-the-end iterator
  OutputIterator output,        //add back-inserter
  double retain_percentage = 5, ///< percentage of points to retain
  double neighbor_radius = 0.25, ///< size of neighbors.
  const unsigned int max_iter_number = 35, ///< number of iterations.
  const bool need_compute_density = true ///< if needed to compute density to   
                                          /// generate more uniform result. 
)
{
  typedef typename value_type_traits<OutputIterator>::type OutputIteratorType;
  return wlop_simplify_and_regularize_point_set
    <Concurrency_tag, OutputIteratorType>(
    first, beyond,
    output,
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    make_dereference_property_map(first),
#else
    make_identity_property_map(typename std::iterator_traits<RandomAccessIterator >::
                               value_type()),
#endif
    retain_percentage, 
    neighbor_radius, 
    max_iter_number, 
    need_compute_density);
}
/// @endcond

} //namespace CGAL

#endif // CGAL_wlop_simplify_and_regularize_point_set_H
