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

#ifndef CGAL_SIMPLIFY_AND_REGULARIZE_POINT_SET_BALL_TREE_H
#define CGAL_SIMPLIFY_AND_REGULARIZE_POINT_SET_BALL_TREE_H

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

#define CGAL_DEBUG_MODE

//#include <tbb/parallel_for.h>
//#include <tbb/blocked_range.h>

namespace CGAL {

/// \cond SKIP_IN_MANUAL

// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace simplify_and_regularize_internal{

// ----------------------------------------------------------------------------
// WLOP algorithm section
// ----------------------------------------------------------------------------

/// Compute anverage and repulsion term and then 
///	compute update sample points positions
/// 
/// \pre `radius > 0`
///
/// @tparam Kernel Geometric traits class.
///
/// @return average term vector
template <typename Kernel>
typename Kernel::Point_3
compute_update_sample_point(
  const typename Kernel::Point_3& query, ///< 3D point to project
  const std::vector<rich_grid_internal::Rich_point<Kernel> >& 
                    neighbor_original_points,///< neighbor original points
  const std::vector<rich_grid_internal::Rich_point<Kernel> >& 
                    neighbor_sample_points, ///< neighbor sample points
  const typename Kernel::FT radius, ///<accept neighborhood radius
  const std::vector<typename Kernel::FT>& original_densities, ///<
  const std::vector<typename Kernel::FT>& sample_densities ///< 
)
{
  CGAL_point_set_processing_precondition(radius > 0);
  CGAL_point_set_processing_precondition(neighbor_original_points.size() >= 1);
  bool is_original_densities_empty = original_densities.empty();
  bool is_sample_densities_empty = sample_densities.empty();

  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;
  typedef rich_grid_internal::Rich_point<Kernel> Rich_point;

  FT radius2 = radius * radius;

  //Compute average term
  FT weight = (FT)0.0, average_weight_sum = (FT)0.0;
  FT iradius16 = -(FT)4.0/radius2;
  Vector average = CGAL::NULL_VECTOR; 

  std::vector<Rich_point>::const_iterator iter;
  iter = neighbor_original_points.begin();
  for (; iter != neighbor_original_points.end(); ++iter)
  {
    const Point& np = iter->pt;
    unsigned int idx_of_original = iter->index;

    FT dist2 = CGAL::squared_distance(query, np);
    weight = exp(dist2 * iradius16);

    if(!is_original_densities_empty)
    {
      weight *= original_densities[idx_of_original];
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
  
  //Compute repulsion term
  weight = (FT)0.0;
  FT repulsion_weight_sum = (FT)0.0;

  Vector repulsion = CGAL::NULL_VECTOR; 
  iter = neighbor_sample_points.begin();
  for (; iter != neighbor_sample_points.end(); ++iter)
  {
	  const Point& np = iter->pt;
	  unsigned int idx_of_sample = iter->index;

	  Vector diff = query - np;

	  FT dist2 = CGAL::squared_distance(query, np);
	  FT dist = std::sqrt(dist2);

	  weight = std::exp(dist2 * iradius16) * std::pow(FT(1.0)/dist, 2);
	  if(!is_sample_densities_empty)
	  {
		  weight *= sample_densities[idx_of_sample];
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

  // Compute update sample point
  Point update_sample = CGAL::ORIGIN + average + FT(0.5) * repulsion;
  return update_sample;
}

/// Compute density weight for each original points,
/// according to their ball neighborhood original points,
/// only need to compute one time in the first iteration
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

  std::vector<Point>::const_iterator iter;
  iter = neighbor_original_points.begin();
  for (; iter != neighbor_original_points.end(); ++iter)
  {
    const Point& np = *iter;
    FT dist2 = CGAL::squared_distance(query, np);
    density_weight += std::exp(dist2 * iradius16);
  }

  // output
  return FT(1.0) / density_weight;
}


/// Compute density weight for sample point,
/// according to their ball neighborhood sample points,
/// need to re-compute in each iterations.
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

  std::vector<Point>::const_iterator iter;
  iter = neighbor_sample_points.begin();
  for (; iter != neighbor_sample_points.end(); ++iter)
  {
    const Point& np = *iter;
    FT dist2 = CGAL::squared_distance(query, np);
    density_weight += std::exp(dist2 * iradius16);
  }

  // output
  return density_weight;
}

/// Extract Points information from indexes, a helper function
///
/// @tparam Kernel Geometric traits class.
///
/// @return Points
template <typename Kernel>
std::vector<typename Kernel::Point_3>
  get_points_from_indexes(
  const std::vector<unsigned int> indexes, ///< indexes
  const std::vector<typename Kernel::Point_3>& all_points ///< all points
  )
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  std::vector<Point> output_points(indexes.size());

  // extract points
  std::vector<Point>::iterator points_iter = output_points.begin();
  unsigned int i = 0;
  for (; points_iter != output_points.end(); ++points_iter)
  {
    *points_iter = all_points[indexes[i++]];
  }

  // output
  return output_points;
}

} // namespace wlop_regularize_and_simplify_point_set

/// \endcond

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
/// More deatail see: http://web.siat.ac.cn/~huihuang/WLOP/WLOP_page.html
///
/// @tparam RandomAccessIterator  iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` 
///         with a value_type = Point_3<Kernel>.
///         It can be omitted if RandomAccessIterator  value_type is convertible to 
///         Point_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///      It can be omitted and deduced automatically from PointPMap value_type.
///
/// @return iterator of the first point to downsampled points.

// This variant requires all parameters.
template <typename RandomAccessIterator , typename PointPMap, typename Kernel>
RandomAccessIterator 
wlop_simplify_and_regularize_point_set(
  RandomAccessIterator  first,  ///< iterator over the first input point.
  RandomAccessIterator  beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map RandomAccessIterator  -> Point_3
  const typename Kernel::FT retain_percentage, ///< percentage to retain,
                                               ///default is 5%
  typename Kernel::FT neighbor_radius, ///< size of neighbors,
                                       /// if the value is negative or non-specific,                                      
                                       /// the estimate value is diameter of bounding box * 0.05.
  const unsigned int max_iter_number,///< number of iterations, rang from 30 to 100 is good,
                                     /// default is 35.
  const bool need_compute_density, ///< if needed to compute density to 
                                   ///generate more rugularized result.                                
  const Kernel& /*kernel*/ ///< geometric traits.
)
{
  Timer task_timer;

  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;
  typedef typename rich_grid_internal::Rich_point<Kernel> Rich_point;

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
  RandomAccessIterator  it;
  RandomAccessIterator  first_original_point = first;
  RandomAccessIterator  first_sample_point = first;
  std::advance(first_sample_point, first_index_to_sample);

  //Copy sample points
  std::vector<Point> sample_points(nb_points_sample);
  std::vector<Point>::iterator sample_iter = sample_points.begin();
  for(it = first_sample_point; it != beyond; ++it, ++sample_iter)
  {
  #ifdef CGAL_USE_PROPERTY_MAPS_API_V1
      *sample_iter = get(point_pmap, it);
  #else
      *sample_iter = get(point_pmap, *it);
  #endif
  }
    
  //Copy original points(Maybe not the best choice)
  std::vector<Point> original_points(nb_points_original);
  std::vector<Point>::iterator original_iter = original_points.begin();
  for(it = first_original_point; it != beyond; ++it, ++original_iter)
  {
  #ifdef CGAL_USE_PROPERTY_MAPS_API_V1
     *original_iter = get(point_pmap, it);
  #else
     *original_iter = get(point_pmap, *it);
  #endif
  }

  // Initialization
  std::vector<Rich_point> original_rich_points(nb_points_original);
  std::vector<Rich_point> sample_rich_points(nb_points_sample);

  CGAL::Bbox_3 bbox(0, 0, 0, 0, 0, 0);
  
  original_iter = original_points.begin();
  int index = 0;
  std::vector<Rich_point>::iterator origianl_rich_iter;
  origianl_rich_iter = original_rich_points.begin();
  for (; original_iter != original_points.end();
         ++original_iter,  ++origianl_rich_iter)
  {  
    *origianl_rich_iter = Rich_point(*original_iter, index++);
    bbox += original_iter->bbox();
  }

  //compute default neighbor_radius
  if (neighbor_radius < 0)
  {
    Point max_p(bbox.xmax(), bbox.ymax(), bbox.zmax());
    Point min_p(bbox.xmin(), bbox.ymin(), bbox.zmin());
    FT bbox_diameter = CGAL::squared_distance(max_p, min_p);
    neighbor_radius = std::sqrt(bbox_diameter) * 0.05;
  
  #ifdef CGAL_DEBUG_MODE
    std::cout << "default estimate radius:  " << neighbor_radius 
              << std::endl << std::endl;
  #endif
  }
  CGAL_point_set_processing_precondition(neighbor_radius > 0);

  // Compute original density weight for original points if user needed
  std::vector<FT> original_densities;
  if (need_compute_density)
  {
  #ifdef CGAL_DEBUG_MODE
    task_timer.start();
    std::cout << "Initialization / Compute Density For Original" << std::endl;
  #endif

    rich_grid_internal::compute_ball_neighbors_one_self(original_rich_points,
                                                        bbox, 
                                                        neighbor_radius);
              
    origianl_rich_iter = original_rich_points.begin();
    original_iter = original_points.begin();
    for (; original_iter != original_points.end(); 
          ++original_iter, ++origianl_rich_iter)
    {
      //get original point positions from indexes
      std::vector<Point> original_neighbors = 
           simplify_and_regularize_internal::get_points_from_indexes<Kernel>
                                            (origianl_rich_iter->neighbors,
                                             original_points );

      FT density = simplify_and_regularize_internal::
                   compute_density_weight_for_original_point<Kernel>(
                   *original_iter,
                   original_neighbors,
                   neighbor_radius);

      original_densities.push_back(density);
      origianl_rich_iter->neighbors.clear();
      //original_rich_points[i].neighbors.swap(std::vector<unsigned int>());
    }

  #ifdef CGAL_DEBUG_MODE
    long memory = CGAL::Memory_sizer().virtual_size();
    std::cout << "done: " << task_timer.time() << " seconds, " 
      << (memory>>20) << " Mb allocated" << std::endl << std::endl;
    task_timer.stop();
  #endif
  }


  // initial rich sample points
  std::vector<Rich_point>::iterator sample_rich_iter;
  sample_rich_iter = sample_rich_points.begin();
  sample_iter = sample_points.begin();
  index = 0;
  for (; sample_iter != sample_points.end(); ++sample_iter, ++sample_rich_iter)
  {
    *sample_rich_iter = Rich_point(*sample_iter, index++);
  }


  std::vector<Point> update_sample_points(nb_points_sample);
  std::vector<Point>::iterator update_sample_iter;
  for (unsigned int iteration = 0; iteration < max_iter_number; iteration++)
  {
    // Build rich-grid for sample-sample neighborhood
    rich_grid_internal::compute_ball_neighbors_one_self(sample_rich_points,
                                                        bbox, 
                                                        neighbor_radius);
  

    // Compute sample density weight for sample points if user needed
    std::vector<FT> sample_densities;
    if (need_compute_density)
    {
	    sample_iter = sample_points.begin();
	    sample_rich_iter = sample_rich_points.begin();
	    for (; sample_rich_iter != sample_rich_points.end(); 
		         ++sample_iter, ++sample_rich_iter)
      {
        std::vector<Point> sample_neighbors = 
        simplify_and_regularize_internal::get_points_from_indexes<Kernel>(
                                          sample_rich_iter->neighbors,
                                          sample_points);
      
        FT density = simplify_and_regularize_internal::
                     compute_density_weight_for_sample_point<Kernel>
                     (*sample_iter, sample_neighbors, neighbor_radius);
      
        sample_densities.push_back(density);
      }
    }

    // Build rich-grid for sample-original neighborhood
    rich_grid_internal::compute_ball_neighbors_one_to_another
                       (sample_rich_points,
                        original_rich_points, bbox, neighbor_radius);
 
	  sample_iter = sample_points.begin();
	  sample_rich_iter = sample_rich_points.begin();
	  for (update_sample_iter = update_sample_points.begin();
	  	   update_sample_iter != update_sample_points.end();
	  	   ++update_sample_iter, ++sample_iter, ++sample_rich_iter)
	  {
	  	Point& p = *sample_iter;
    
	  	std::vector<Rich_point> rich_original_neighbors;
	  	std::vector<unsigned int>& original_neighors_indexes = 
	  		                       sample_rich_iter->original_neighbors;
    
      std::vector<unsigned int>::iterator iter;
	  	for (iter = original_neighors_indexes.begin();
           iter != original_neighors_indexes.end(); ++iter)
	  	{
	  		unsigned int idx_of_original = *iter;
	  		Rich_point rp(original_points[idx_of_original], idx_of_original);
	  		rich_original_neighbors.push_back(rp);
	  	}
    
	  	std::vector<Rich_point> rich_sample_neighbors;
	  	std::vector<unsigned int>& sample_neighors_indexes = 
	  		                         sample_rich_iter->neighbors;
    
      for (iter = sample_neighors_indexes.begin();
           iter != sample_neighors_indexes.end(); ++iter)
	  	{
	  	  unsigned int idx_of_sample = *iter;
	  	  Rich_point rp(sample_points[idx_of_sample], idx_of_sample);
	  	  rich_sample_neighbors.push_back(rp);
	  	}
    
	  	*update_sample_iter = simplify_and_regularize_internal::
                            compute_update_sample_point<Kernel>
                            (*sample_iter,
	  						         	   rich_original_neighbors,
	  							           rich_sample_neighbors,
                             neighbor_radius,
                             original_densities,
                             sample_densities);
	  }


    // update sample points positions
    sample_iter = sample_points.begin();
    sample_rich_iter = sample_rich_points.begin();
    update_sample_iter = update_sample_points.begin();
    for (; sample_iter != sample_points.end(); 
           ++sample_iter, ++sample_rich_iter, ++update_sample_iter)
    {
      *sample_iter = *update_sample_iter;
      sample_rich_iter->pt = *sample_iter;
    }
  
  #ifdef CGAL_DEBUG_MODE
    std::cout << "Compute average & repulsion term and updated" << std::endl;
    long memory = CGAL::Memory_sizer().virtual_size();
    std::cout << "done: " << task_timer.time() << " seconds, " 
      << (memory>>20) << " Mb allocated" << std::endl;
    task_timer.stop();
    task_timer.start();
    std::cout << "iterate:  " << iteration + 1 << "	"<< std::endl << std::endl;
  #endif
  }

  //Copy back modified sample points to original points for output
  sample_iter = sample_points.begin();
  for(it = first_sample_point; it != beyond; ++it, ++sample_iter)
  {
    Point& sample_p = *sample_iter;

  #ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    Point& original_p = get(point_pmap, it);
    original_p = sample_p;
  #else
    Point& original_p = get(point_pmap, *it);
    original_p = sample_p;
  #endif
  }

  #ifdef CGAL_DEBUG_MODE
    original_rich_points.erase(original_rich_points.begin(), 
                                  original_rich_points.end());
    
    original_rich_points.clear();

    std::cout << "Copy back done: " << task_timer.time() 
      << " seconds "  << std::endl;

    task_timer.stop();
    task_timer.start();
    original_rich_points.swap(std::vector<Rich_point>());
    
    std::cout << "STL release memory: " << task_timer.time() 
              << " seconds "  << std::endl;
    task_timer.stop();
    
    sample_rich_points.clear();
    sample_points.clear();
  #endif

  return first_sample_point;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename RandomAccessIterator , typename PointPMap>
RandomAccessIterator 
wlop_simplify_and_regularize_point_set(
  RandomAccessIterator  first, ///< iterator over the first input point
  RandomAccessIterator  beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map RandomAccessIterator  -> Point_3
  double retain_percentage, ///< percentage of points to retain
  double neighbor_radius, ///< size of neighbors.
  const unsigned int max_iter_number, ///< number of iterations.
  const bool need_compute_density  ///< if needed to compute density 
                                   ///  to generate more rugularized result.                                 
) 
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return wlop_simplify_and_regularize_point_set(
    first, beyond,
    point_pmap,
    retain_percentage,
    neighbor_radius,
    max_iter_number,
    need_compute_density,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
/// This variant creates a default point property map=Dereference_property_map.
template <typename RandomAccessIterator >
RandomAccessIterator 
wlop_simplify_and_regularize_point_set(
  RandomAccessIterator  first, ///< iterator over the first input point
  RandomAccessIterator  beyond, ///< past-the-end iterator
  double retain_percentage = 5, ///< percentage of points to retain
  double neighbor_radius = -1, ///< size of neighbors.
  const unsigned int max_iter_number = 35, ///< number of iterations.
  const bool need_compute_density = true ///< if needed to compute density to   
                                          /// generate more uniform result. 
)
{
  return wlop_simplify_and_regularize_point_set(
    first, beyond,
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

#endif // CGAL_SIMPLIFY_AND_REGULARIZE_POINT_SET_BALL_TREE_H
