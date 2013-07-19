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
// Author(s) : Shihao Wu, Cl¨¦ment Jamin, Pierre Alliez 

#ifndef CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_BALL_TREE_H
#define CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_BALL_TREE_H

#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Rich_grid.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>
#include <iterator>

#include <boost/version.hpp>
#if BOOST_VERSION >= 104000
  #include <boost/property_map/property_map.hpp>
#else
  #include <boost/property_map.hpp>
#endif

//#include <tbb/parallel_for.h>
//#include <tbb/blocked_range.h>

namespace CGAL {

/// \cond SKIP_IN_MANUAL

// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace regularize_and_simplify_internal{



// ----------------------------------------------------------------------------
// WLOP algorithm section
// ----------------------------------------------------------------------------
/// Compute average term for each sample points
/// According to their ball neighborhood original points
/// 
/// \pre `radius > 0`
///
/// @tparam Kernel Geometric traits class.
///
/// @return average term vector
template <typename Kernel>
typename Kernel::Vector_3
compute_average_term(
  const typename Kernel::Point_3& query, ///< 3D point to project
  const std::vector<rich_grid_internel::Rich_point<Kernel> >& 
                    neighbor_original_points,///< neighbor sample points
  const typename Kernel::FT radius, ///<accept neighborhood radius
  const std::vector<typename Kernel::FT>& density_weight_set ///< densities
)
{
  CGAL_point_set_processing_precondition(radius > 0);
  CGAL_point_set_processing_precondition(neighbor_original_points.size() >= 1);
  bool is_density_weight_set_empty = density_weight_set.empty();

  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;

  FT radius2 = radius * radius;

  //Compute average term
  FT weight = (FT)0.0, average_weight_sum = (FT)0.0;
  FT iradius16 = -(FT)4.0/radius2;
  Vector average = CGAL::NULL_VECTOR; 
  for (unsigned int i = 0; i < neighbor_original_points.size(); i++)
  {
    const Point& np = neighbor_original_points[i].pt;
    unsigned int idx_of_original = neighbor_original_points[i].index;

    FT dist2 = CGAL::squared_distance(query, np);
    weight = exp(dist2 * iradius16);

    if(!is_density_weight_set_empty)
    {
      weight *= density_weight_set[idx_of_original];
    }

    average_weight_sum += weight;
    average = average + (np - CGAL::ORIGIN) * weight;
  }

  // output
  return average/average_weight_sum;
}

/// Compute repulsion term for each sample points
/// According to their Ball neighborhood sample points
/// 
/// \pre \pre `radius > 0`
///
/// @tparam Kernel Geometric traits class.
///
/// @return repulsion term vector
template <typename Kernel>
typename Kernel::Vector_3
compute_repulsion_term(
  const typename Kernel::Point_3& query, ///< 3D point to project
  const std::vector<rich_grid_internel::Rich_point<Kernel> >& 
             neighbor_sample_points, ///< neighbor sample points
  const typename Kernel::FT radius, ///<accept neighborhood radius
  const std::vector<typename Kernel::FT>& density_weight_set ///< densities
)
{
  CGAL_point_set_processing_precondition(radius > 0);
  CGAL_point_set_processing_precondition(neighbor_sample_points.size() >= 1);

  bool is_density_weight_set_empty = density_weight_set.empty();

  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;

  FT radius2 = radius * radius;

  //Compute average term
  FT weight = (FT)0.0, repulsion_weight_sum = (FT)0.0;
  FT iradius16 = -(FT)4.0/radius2;

  Vector repulsion = CGAL::NULL_VECTOR; 
  for (unsigned int i = 0; i < neighbor_sample_points.size(); i++)
  {
    const Point& np = neighbor_sample_points[i].pt;
    unsigned int idx_of_sample = neighbor_sample_points[i].index;

    Vector diff = query - np;

    FT dist2 = CGAL::squared_distance(query, np);
    FT dist = std::sqrt(dist2);

    weight = std::exp(dist2 * iradius16) * std::pow(FT(1.0)/dist, 2);
    if(!is_density_weight_set_empty)
    {
      weight *= density_weight_set[idx_of_sample];
    }

    repulsion_weight_sum += weight;
    repulsion = repulsion + diff * weight;
  }

  // output
  return repulsion/repulsion_weight_sum;
}




/// Compute density weight for each original points,
/// according to their ball neighborhood original points
/// 
/// \pre `radius > 0`
///
/// @tparam Kernel Geometric traits class.
///
/// @return density weight
template <typename Kernel>
typename Kernel::FT
compute_density_weight_for_original_point(
  const typename Kernel::Point_3& query, ///< 3D point to project
  const std::vector<typename Kernel::Point_3>& neighbor_original_points, ///< 
  const typename Kernel::FT radius ///<accept neighborhood radius
)
{
  CGAL_point_set_processing_precondition(radius > 0);

  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;

  //Compute density weight
  FT radius2 = radius * radius;
  FT density_weight = (FT)1.0;
  FT iradius16 = -(FT)4.0/radius2;

  for (unsigned int i = 0; i < neighbor_original_points.size(); i++)
  {
    const Point& np = neighbor_original_points[i];
    FT dist2 = CGAL::squared_distance(query, np);
    density_weight += std::exp(dist2 * iradius16);
  }

  // output
  return FT(1.0) / density_weight;
}


/// Compute density weight for sample point,
/// according to their ball neighborhood sample points
/// 
/// \pre `radius > 0`
///
/// @tparam Kernel Geometric traits class.
///
/// @return density weight
template <typename Kernel>
typename Kernel::FT
compute_density_weight_for_sample_point(
  const typename Kernel::Point_3& query, ///< 3D point to project
  const std::vector<typename Kernel::Point_3>& neighbor_sample_points, ///< 
  const typename Kernel::FT radius
)
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;

  //Compute density weight
  FT radius2 = radius * radius;
  FT density_weight = (FT)1.0;
  FT iradius16 = -(FT)4.0/radius2;

  for (unsigned int i = 0; i < neighbor_sample_points.size(); i++)
  {
    const Point& np = neighbor_sample_points[i];
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
///	denoised, outlier-free and evenly distributed particles over the original 
/// dense point cloud, so as to improve the reliability of other algorithms. 
///
///	The core of the algorithm is a Weighted Locally Optimal projection operator
/// with a density uniformization term. 
/// More deatail please see:http://web.siat.ac.cn/~huihuang/WLOP/WLOP_page.html
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
wlop_regularize_and_simplify_point_set_using_rich_grid(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  const typename Kernel::FT retain_percentage, ///< percentage to retain.
  const typename Kernel::FT neighbor_radius, ///< size of neighbors.
  const unsigned int iter_number,///< number of iterations.
  const bool need_compute_density, ///< if needed to compute density to 
                                   ///generate more rugularized result.                                
  const Kernel& /*kernel*/ ///< geometric traits.
)
{
  CGAL_point_set_processing_precondition(neighbor_radius > 0);
  Timer task_timer;

  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;
  typedef typename rich_grid_internel::Rich_point<Kernel> Rich_point;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(first != beyond);
  CGAL_point_set_processing_precondition(retain_percentage >= 0 && 
                                         retain_percentage <= 100);

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
  for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
  {
  #ifdef CGAL_USE_PROPERTY_MAPS_API_V1
      sample_points[i] = get(point_pmap, *it);
  #else
      sample_points[i] = get(point_pmap, *it);
  #endif
  }
    
  //Copy original points(Maybe not the best choice)
  std::vector<Point> original_points(nb_points_original);
  for(it = first_original_point, i = 0; it != beyond; ++it, i++)
  {
  #ifdef CGAL_USE_PROPERTY_MAPS_API_V1
     original_points[i] = get(point_pmap, it);
  #else
     original_points[i] = get(point_pmap, *it);
  #endif
  }

  // Initilization
  std::vector<Rich_point> original_rich_point_set(nb_points_original);
  std::vector<Rich_point> sample_rich_point_set(nb_points_sample);

  rich_grid_internel::Rich_box<Kernel> box;
  for (i = 0; i < nb_points_original; i++)
  {  
    Point& p0 = original_points[i];
    Rich_point rp(p0, i);
    original_rich_point_set[i] = rp;
    box.add_point(rp.pt);
  }

  // Compute original density weight for original points if user needed
  std::vector<FT> original_density_weight_set;
  if (need_compute_density)
  {
    task_timer.start();
    std::cout << "Initialization / Compute Density For Original" << std::endl;
  
    rich_grid_internel::compute_ball_neighbors_one_self(original_rich_point_set
                                                       , box, neighbor_radius);
                        
    for (i = 0; i < nb_points_original; i++)
    {
      // get original point positions from indexes
      std::vector<Point> original_neighbors;
      std::vector<unsigned int>& neighors_indexes = 
                                 original_rich_point_set[i].neighbors;

      for (unsigned int j = 0; j < neighors_indexes.size(); j++)
      {
        original_neighbors.push_back(original_points[neighors_indexes[j]]);
      }

      // compute density
      FT density = regularize_and_simplify_internal::
                   compute_density_weight_for_original_point<Kernel>(
                   original_points[i],
                   original_neighbors, neighbor_radius);

      original_density_weight_set.push_back(density);
      original_rich_point_set[i].neighbors.clear();
    }

    long memory = CGAL::Memory_sizer().virtual_size();
    std::cout << "done: " << task_timer.time() << " seconds, " 
      << (memory>>20) << " Mb allocated" << std::endl << std::endl;
    task_timer.stop();
  }


  for (unsigned int iter_n = 0; iter_n < iter_number; iter_n++)
  {
    task_timer.start();
    std::cout << "Compute average term and repulsion term " << std::endl;

    // Build Ball Tree For Sample Neighbor
    for (i=0 ; i < sample_points.size(); i++)
    {
      Point& p0 = sample_points[i];
      Rich_point rp(p0, i);
      sample_rich_point_set[i] = rp;
    }
    rich_grid_internel::compute_ball_neighbors_one_self(sample_rich_point_set,
                                                        box, neighbor_radius);

    // Compute sample density weight for sample points if user needed
    std::vector<FT> sample_density_weight_set;
    if (need_compute_density)
    {
      for (i=0 ; i < sample_points.size(); i++)
      {
        std::vector<Point> sample_neighbors;
        std::vector<unsigned int>& neighors_indexes = 
                                   sample_rich_point_set[i].neighbors;

        for (unsigned int j = 0; j < neighors_indexes.size(); j++)
        {
          sample_neighbors.push_back(sample_points[neighors_indexes[j]]);
        }

        FT density = regularize_and_simplify_internal::
                     compute_density_weight_for_sample_point<Kernel>
                     (sample_points[i], sample_neighbors, neighbor_radius);

        sample_density_weight_set.push_back(density);
      }
    }

    // Build Ball Tree For Sample-Original Neighbor
    rich_grid_internel::compute_ball_neighbors_one_to_another
                       (sample_rich_point_set,
                        original_rich_point_set, box, neighbor_radius);

    // Compute average term and repulsion term for each sample points,
    // then update each sample points
    std::vector<Vector> average_set(nb_points_sample);
    std::vector<Vector> repulsion_set(nb_points_sample);

    // average term
    for (i = 0; i < sample_points.size(); i++)
    {
      Point& p = sample_points[i];
      std::vector<Rich_point> rich_original_neighbors;
      std::vector<unsigned int>& neighors_indexes = 
                                 sample_rich_point_set[i].original_neighbors;

      if (neighors_indexes.empty())
      {
        average_set[i] = p - CGAL::ORIGIN;
        continue;
      }

      for (unsigned int j = 0; j < neighors_indexes.size(); j++)
      {
        unsigned int idx_of_original = neighors_indexes[j];
        Rich_point rp(original_points[idx_of_original], idx_of_original);
        rich_original_neighbors.push_back(rp);
      }

      average_set[i] = regularize_and_simplify_internal::
                       compute_average_term<Kernel>
                       (p, rich_original_neighbors, neighbor_radius, 
                        original_density_weight_set); 
    }

    //repulsion term
    for (i = 0; i < sample_points.size(); i++)
    {
      std::vector<Rich_point> rich_sample_neighbors;
      std::vector<unsigned int>& neighors_indexes = 
                                 sample_rich_point_set[i].neighbors;
        
      if (neighors_indexes.empty())
      {
        repulsion_set[i] = CGAL::NULL_VECTOR;
        continue;
      }

      for (unsigned int j = 0; j < neighors_indexes.size(); j++)
      {
        unsigned int idx_of_sample = neighors_indexes[j];
        Rich_point rp(sample_points[idx_of_sample], idx_of_sample);
        rich_sample_neighbors.push_back(rp);
      }

      Point& p = sample_points[i];
      repulsion_set[i] = regularize_and_simplify_internal::
                         compute_repulsion_term<Kernel>
                         (p, rich_sample_neighbors, neighbor_radius, 
                         sample_density_weight_set);
    }
    
    // update points positions according to average and repulsion term
    for (i = 0; i < sample_points.size(); i++)
    {
      Point& p = sample_points[i];
      p = CGAL::ORIGIN + average_set[i] + (FT)0.5 * repulsion_set[i];
    }

    long memory = CGAL::Memory_sizer().virtual_size();
    std::cout << "done: " << task_timer.time() << " seconds, " 
      << (memory>>20) << " Mb allocated" << std::endl;
    task_timer.stop();

    std::cout << "iterate:	" << iter_n + 1 <<  "	"<< std::endl << std::endl;
  }

  //Copy back modified sample points to original points for output
  for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
  {
    Point& sample_p = sample_points[i];

  #ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    Point& original_p = get(point_pmap, it);
    original_p = sample_p;
  #else
    Point& original_p = get(point_pmap, *it);
    original_p = sample_p;
  #endif

  //#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
  //    put(point_pmap, sample_p, it);
  //#else
  //    put(point_pmap, sample_p, *it);
  //#endif
  }

  return first_sample_point;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename ForwardIterator, typename PointPMap>
ForwardIterator
wlop_regularize_and_simplify_point_set_using_rich_grid(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double retain_percentage, ///< percentage of points to retain
  double neighbor_radius, ///< size of neighbors.
  const unsigned int iter_number, ///< number of iterations.
  const bool need_compute_density  ///< if needed to compute density 
                                   ///  to generate more rugularized result.                                 
) 
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return wlop_regularize_and_simplify_point_set_using_rich_grid(
    first, beyond,
    point_pmap,
    retain_percentage,
    neighbor_radius,
    iter_number,
    need_compute_density,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Dereference_property_map.
template <typename ForwardIterator>
ForwardIterator
wlop_regularize_and_simplify_point_set_using_rich_grid(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  double retain_percentage, ///< percentage of points to retain
  double neighbor_radius, ///< size of neighbors.
  const unsigned int iter_number, ///< number of iterations.
  const bool need_compute_density ///< if needed to compute density to generate
                                  ///  more rugularized result. 
)
{
  return wlop_regularize_and_simplify_point_set_using_rich_grid(
    first, beyond,
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    make_dereference_property_map(first),
#else
    make_identity_property_map(typename std::iterator_traits<ForwardIterator>::
                               value_type()),
#endif
    retain_percentage, neighbor_radius, iter_number, need_compute_density);
}
/// @endcond

} //namespace CGAL

#endif // CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H
