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

#ifndef CGAL_UPSAMPLE_POINT_SET_H
#define CGAL_UPSAMPLE_POINT_SET_H

#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/internal/Rich_grid.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/compute_average_spacing.h>

#include <iterator>
#include <set>
#include <utility>

#include <boost/version.hpp>
#if BOOST_VERSION >= 104000
  #include <boost/property_map/property_map.hpp>
#else
  #include <boost/property_map.hpp>
#endif

//#define  CGAL_PSP3_VERBOSE

namespace CGAL {

// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL

namespace upsample_internal{

/// For each query point, select a best "base point" in its neighborhoods.
/// Then, a new point will be interpolated between query point and "base point".
/// This is the key part of the upsample algorithm 
/// 
/// \pre `radius > 0`
///
/// @tparam Kernel Geometric traits class.
///
/// @return local density length
template <typename Kernel>
typename Kernel::FT
base_point_selection(
  const rich_grid_internal::Rich_point<Kernel>& query, ///< 3D point to project
  const std::vector<rich_grid_internal::Rich_point<Kernel> >& 
                    neighbor_points,///< neighbor sample points
  const typename Kernel::FT edge_sensitivity,///< edge senstivity parameter
  unsigned int& output_base_index ///< base point index
  )
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;
  typedef typename rich_grid_internal::Rich_point<Kernel> Rich_point;

  if (neighbor_points.empty())
  {
#ifdef CGAL_PSP3_VERBOSE
    std::cout << "empty neighborhood" << std::endl;
#endif
    output_base_index = query.index;
    return 0.0;
  }

  FT best_dist2 = -10.0;
  const Rich_point& v = query;
  typename std::vector<Rich_point>::const_iterator iter = neighbor_points.begin();
  for (; iter != neighbor_points.end(); ++iter)
  {
    const Point& t = iter->pt;

    const Vector& vm = v.normal;
    const Vector& tm = iter->normal;

    Vector diff_v_t = t - v.pt;
    Point mid_point = v.pt + (diff_v_t * FT(0.5));
    
    FT dot_produce = std::pow((FT(2.0) - vm * tm), edge_sensitivity);

    Vector diff_t_mid = mid_point - t;
    FT project_t = diff_t_mid * tm;
    FT min_dist2 = diff_t_mid.squared_length() - project_t * project_t;

    typename std::vector<Rich_point>::const_iterator iter_in = neighbor_points.begin();
    for (; iter_in != neighbor_points.end(); ++iter_in)
    {
      Vector diff_s_mid = mid_point - iter_in->pt;
      FT project_s = diff_s_mid * iter_in->normal;

      FT proj_min2 = diff_s_mid.squared_length() - project_s * project_s;

      if (proj_min2 < min_dist2)
      {
        min_dist2 = proj_min2;
      }
    }
    min_dist2 *= dot_produce;

    if (min_dist2 > best_dist2)
    {
      best_dist2 = min_dist2;
      output_base_index = iter->index;
    }
  }

  return best_dist2; 
}

/// For each new inserted point, we need to do the following job
/// 1, get neighbor information from the two "parent points"
/// 2, update position and determine normal by bilateral projection 
/// 3, update neighbor information again
///
/// \pre `radius > 0`
///
/// @tparam Kernel Geometric traits class.
///
template <typename Kernel>
void
update_new_point(
  unsigned int new_point_index, ///< new inserted point
  unsigned int father_index, ///< father point index
  unsigned int mother_index, ///< mother point index
  std::vector<rich_grid_internal::Rich_point<Kernel> >& rich_point_set,
                                                           ///< all rich points
  const typename Kernel::FT radius,          ///< accept neighborhood radius
  const typename Kernel::FT sharpness_bandwidth  ///< control sharpness
)
{
  // basic geometric types
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;
  typedef typename rich_grid_internal::Rich_point<Kernel> Rich_point;

  CGAL_assertion_code( unsigned int size = static_cast<unsigned int>(rich_point_set.size()) );
  CGAL_point_set_processing_precondition(father_index >= 0 &&
                                         father_index < size);
  CGAL_point_set_processing_precondition(mother_index >= 0 &&
                                           mother_index < size);

  // 1, get neighbor information from the two "parent points"
  Rich_point& new_v = rich_point_set[new_point_index];
  Rich_point& father_v = rich_point_set[father_index];
  Rich_point& mother_v = rich_point_set[mother_index];

  std::set<int> neighbor_indexes;
  std::vector<unsigned int>::iterator iter;
  for (iter = father_v.neighbors.begin();
       iter != father_v.neighbors.end();
       ++iter)
  {
    neighbor_indexes.insert(*iter);
  }

  for (iter = mother_v.neighbors.begin();
       iter != mother_v.neighbors.end();
       ++iter)
  {
    neighbor_indexes.insert(*iter);
  }

  neighbor_indexes.insert(father_v.index);
  neighbor_indexes.insert(mother_v.index);

  double radius2 = radius * radius;

  new_v.neighbors.clear();
  std::set<int>::iterator set_iter;
  for (set_iter = neighbor_indexes.begin(); 
       set_iter != neighbor_indexes.end(); ++set_iter)
  {
    Rich_point& t = rich_point_set[*set_iter];
    FT dist2 =  CGAL::squared_distance(new_v.pt, t.pt);

    if (dist2 < radius2)
    {
      new_v.neighbors.push_back(t.index);
    }
  }

  // 2, update position and normal by bilateral projection 
  const unsigned int candidate_num = 2; // we have two normal candidates:
                                        // we say father's is 0
                                        //        mother's is 1
  std::vector<Vector> normal_cadidate(candidate_num);
  normal_cadidate[0] = father_v.normal;
  normal_cadidate[1] = mother_v.normal;

  std::vector<FT> project_dist_sum(candidate_num, FT(0.0));
  std::vector<FT> weight_sum(candidate_num, FT(0.0));
  std::vector<Vector> normal_sum(candidate_num, NULL_VECTOR);
   
  FT radius16 = FT(-4.0) / radius2;

  for (unsigned int i = 0; i < new_v.neighbors.size(); ++i)
  {
    Rich_point& t = rich_point_set[new_v.neighbors[i]];
    FT dist2 = CGAL::squared_distance(new_v.pt, t.pt);
    FT theta = std::exp(dist2 * radius16);

    for (unsigned int j = 0; j < candidate_num; j++)
    {
      FT psi = std::exp(-std::pow(1 - normal_cadidate[j] * t.normal, 2)
                       / sharpness_bandwidth);
      FT project_diff_t_v = (t.pt - new_v.pt) * t.normal;
      FT weight = psi * theta;

      project_dist_sum[j] += project_diff_t_v * weight;
      normal_sum[j] = normal_sum[j] + t.normal * weight;
      weight_sum[j] += weight;
    }
  }

  // select best candidate
  FT min_project_dist = (FT)(std::numeric_limits<double>::max)();
  unsigned int best = 0;

  for (unsigned int i = 0; i < candidate_num; ++i)
  {
    FT absolute_dist = abs(project_dist_sum[i] / weight_sum[i]);
    if (absolute_dist < min_project_dist)
    {
      min_project_dist = absolute_dist;
      best = i;
    }
  }

  // update position and normal
  Vector update_normal = normal_sum[best] / weight_sum[best];
  new_v.normal = update_normal / sqrt(update_normal.squared_length());

  FT project_dist = project_dist_sum[best] / weight_sum[best];
  new_v.pt = new_v.pt + new_v.normal * project_dist;


  // 3, update neighbor information again
  new_v.neighbors.clear();
  for (set_iter = neighbor_indexes.begin(); 
       set_iter != neighbor_indexes.end(); ++set_iter)
  {
    Rich_point& t = rich_point_set[*set_iter];
    FT dist2 =  CGAL::squared_distance(new_v.pt, t.pt);

    if (dist2 < radius2)
    {
      new_v.neighbors.push_back(t.index);
      t.neighbors.push_back(new_v.index);
    }
  }
}

} /* namespace upsample_internal */

/// \endcond

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/// \ingroup PkgPointSetProcessing
/// This method progressively upsamples the point set while 
/// approaching the edge singularities (detected by normal variation), which 
/// generates a denser point set from an input point set. This has applications 
/// in point-based rendering, hole filling, and sparse surface reconstruction. 
/// Normals of points are required as input. For more details, please refer to \cgalCite{ear-2013}.
/// 
/// @tparam OutputIterator Type of the output iterator. 
///         The type of the objects put in it is 
///         `std::pair<Kernel::Point_3, Kernel::Vector_3>`.
///         Note that the user may use a 
///         <A HREF="http://www.boost.org/libs/iterator/doc/function_output_iterator.html">function_output_iterator</A>
///         to match specific needs.
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` 
///         with the value type of `ForwardIterator` as key and `Kernel::Point_3` as value type.
///         It can be omitted if the value type of `ForwardIterator` is convertible to 
///         `Kernel::Point_3`.
/// @tparam NormalPMap is a model of `ReadablePropertyMap` with the value type of `ForwardIterator` as key
///         and `Kernel::Vector_3` as value type.
/// @tparam Kernel Geometric traits class.
///      It can be omitted and deduced automatically from the value type of  `PointPMap`
///      using `Kernel_traits`.
///

// This variant requires all parameters.
template <typename OutputIterator,
          typename ForwardIterator, 
          typename PointPMap, 
          typename NormalPMap,
          typename Kernel>
OutputIterator
edge_aware_upsample_point_set(
  ForwardIterator first,  ///< forward iterator on the first input point.
  ForwardIterator beyond, ///< past-the-end iterator.
  OutputIterator output,  ///< output iterator where output points and normals 
                          ///< are put.
  PointPMap point_pmap,   ///< point property map.
  NormalPMap normal_pmap, ///< vector property map.
  const typename Kernel::FT sharpness_angle,  ///< 
                    ///< controls the preservation of sharp features. 
                    ///< The larger the value is,
                    ///< the smoother the result will be.
                    ///< The range of possible value is [0, 90].
                    ///< See section \ref Point_set_processing_3Upsample_Parameter2
                    ///< for an example.
  typename Kernel::FT edge_sensitivity,  ///<  
                    ///< larger values of edge-sensitivity give higher priority 
                    ///< to inserting points along sharp features.
                    ///< The range of possible values is `[0, 1]`.
                    ///< See section \ref Point_set_processing_3Upsample_Parameter1
                    ///< for an example.
  typename Kernel::FT neighbor_radius, ///< 
                    ///< indicates the radius of the largest hole that should be filled.
                    ///< The default value is set to 3 times the average spacing of the point set.
                    ///< If the value given by user is smaller than the average spacing, 
                    ///< the function will use the default value instead.
  const unsigned int number_of_output_points,///< number of output
                                             ///< points to generate.
  const Kernel& /*kernel*/ ///< geometric traits.
)
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;
  typedef typename rich_grid_internal::Rich_point<Kernel> Rich_point;

  // preconditions
  CGAL_point_set_processing_precondition(first != beyond);
  CGAL_point_set_processing_precondition(sharpness_angle >= 0 
                                       &&sharpness_angle <= 90);
  CGAL_point_set_processing_precondition(edge_sensitivity >= 0 
                                       &&edge_sensitivity <= 1);
  CGAL_point_set_processing_precondition(neighbor_radius > 0);

  edge_sensitivity *= 10;  // just project [0, 1] to [0, 10].

  std::size_t number_of_input = std::distance(first, beyond);
  CGAL_point_set_processing_precondition(number_of_output_points > number_of_input);


  const unsigned int nb_neighbors = 6; // 1 ring
  FT average_spacing = CGAL::compute_average_spacing(
                   first, beyond,
                   point_pmap,
                   nb_neighbors);

  if (neighbor_radius < average_spacing * 1.0)
  {
    neighbor_radius = average_spacing * 3.0;
#ifdef CGAL_PSP3_VERBOSE
    std::cout << "neighbor radius: " << neighbor_radius << std::endl;
#endif
  }
  
  Timer task_timer;

  // copy rich point set
  std::vector<Rich_point> rich_point_set(number_of_input);
  CGAL::Bbox_3 bbox(0., 0., 0., 0., 0., 0.);
  
  ForwardIterator it = first; // point iterator
  for(unsigned int i = 0; it != beyond; ++it, ++i)
  {
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    rich_point_set[i].pt = get(point_pmap, it);
    rich_point_set[i].normal = get(normal_pmap, it);
#else
    rich_point_set[i].pt = get(point_pmap, *it);
    rich_point_set[i].normal = get(normal_pmap, *it);
#endif

    rich_point_set[i].index = i;
    bbox += rich_point_set[i].pt.bbox();
    CGAL_point_set_processing_precondition(rich_point_set[i].normal.squared_length() > 1e-10);
  }

  // compute neighborhood
  rich_grid_internal::compute_ball_neighbors_one_self(rich_point_set,
                                                      bbox,
                                                      neighbor_radius);

  //
  FT cos_sigma = std::cos(sharpness_angle / 180.0 * 3.1415926);
  FT sharpness_bandwidth = std::pow((CGAL::max)((FT)1e-8, (FT)1.0 - cos_sigma), 2);

  FT sum_density = 0.0;
  unsigned int count_density = 1;
  double max_iter_time = 20;
  FT current_radius = neighbor_radius;
  FT density_pass_threshold = 0.0;

  for (unsigned int iter_time = 0; iter_time < max_iter_time; ++iter_time)
  {
  #ifdef CGAL_PSP3_VERBOSE
     std::cout << std::endl << "iter_time: " << iter_time + 1  << std::endl;
  #endif
    if (iter_time > 0)
    {
      current_radius *= 0.75;
      if (current_radius < density_pass_threshold * 3) //3 ring
      {
        current_radius = density_pass_threshold * 3;
      }
      rich_grid_internal::compute_ball_neighbors_one_self(rich_point_set,
                                                          bbox,
                                                          current_radius);
    }
 #ifdef CGAL_PSP3_VERBOSE
    std::cout << "current radius: " << current_radius << std::endl; 
 #endif

    std::size_t current_size = rich_point_set.size();
    std::vector<bool> is_pass_threshold(current_size, false);

    if (iter_time == 0)
    {
      //estimate density threshold for the first time
      for (unsigned int i = 0; i < rich_point_set.size() * 0.05; ++i)
      {
        Rich_point& v = rich_point_set[i];

        if (v.neighbors.empty())
          continue;

        // extract neighbor rich points by index
        std::vector<Rich_point> neighbor_rich_points(v.neighbors.size());
        for (unsigned int n = 0; n < v.neighbors.size(); n++)
        {
          neighbor_rich_points[n] = rich_point_set[v.neighbors[n]];
        }

        unsigned int base_index = 0;
        double density2 = upsample_internal::
                              base_point_selection(v,
                                                   neighbor_rich_points,
                                                   edge_sensitivity,
                                                   base_index);

        if (density2 < 0)
        {
          continue;
        }

        sum_density += density2;
        count_density++;
      }
    }

    density_pass_threshold = sqrt(sum_density / count_density) * 0.65;
    sum_density = 0.;
    count_density = 1;

    FT density_pass_threshold2 = density_pass_threshold * 
                                 density_pass_threshold;
 #ifdef CGAL_PSP3_VERBOSE
    std::cout << "pass_threshold:  " << density_pass_threshold << std::endl;
 #endif
    // insert new points until all the points' density pass the threshold
    unsigned int max_loop_time = 3;
    unsigned int loop = 0;
    while (true)
    {
   #ifdef CGAL_PSP3_VERBOSE
      std::cout << "loop_time: " << loop + 1 << std::endl;
   #endif
      unsigned int count_not_pass = 0;
      loop++;
      for (unsigned int i = 0; i < rich_point_set.size(); ++i)
      {
        if (is_pass_threshold[i])
        {
          continue;
        }

        Rich_point& v = rich_point_set[i];

        if (v.neighbors.empty())
          continue;

        // extract neighbor rich points by index
        std::vector<Rich_point> neighbor_rich_points(v.neighbors.size());
        for (unsigned int n = 0; n < v.neighbors.size(); ++n)
        {
          neighbor_rich_points[n] = rich_point_set[v.neighbors[n]];
        }

        // select base point 
        unsigned int base_index = 0;
        FT density2 = upsample_internal::
                              base_point_selection(v,
                                                   neighbor_rich_points,
                                                   edge_sensitivity,
                                                   base_index);

        // test if it pass the density threshold
        if (density2 < density_pass_threshold2)
        {
          is_pass_threshold[i] = true;
          continue;
        }
        count_not_pass++;

        sum_density += density2;
        count_density++;

        // insert a new rich point
        unsigned int father_index = v.index;
        unsigned int mother_index = base_index;

        Rich_point new_v;
        Rich_point& base = rich_point_set[mother_index];
        Vector diff_v_base = base.pt - v.pt;
        new_v.pt = v.pt + (diff_v_base * FT(0.5));
        new_v.index = static_cast<unsigned int>(rich_point_set.size());

        unsigned int new_point_index = new_v.index;
        rich_point_set.push_back(new_v);
        is_pass_threshold.push_back(false);

        //update new rich point
        upsample_internal::update_new_point(new_point_index, 
                                            father_index, 
                                            mother_index, 
                                            rich_point_set, 
                                            current_radius,
                                            sharpness_bandwidth);

        if (rich_point_set.size() >= number_of_output_points)
        {
          break;
        }
      }
   #ifdef CGAL_PSP3_VERBOSE
      std::cout << "current size: " << rich_point_set.size() << std::endl;
   #endif
      if (count_not_pass == 0 || 
          loop >= max_loop_time || 
          rich_point_set.size() >= number_of_output_points)
      {
        break;
      }

    }

    if (rich_point_set.size() >= number_of_output_points)
    {
      break;
    }
  }

  for (std::size_t i = number_of_input; i < rich_point_set.size(); ++i)
  {
    Rich_point& v = rich_point_set[i];
    Point point = v.pt;
    Vector normal = v.normal;
    *output++ = std::make_pair(point, normal);
  }
 
  return output;
}


/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
template <typename OutputIterator,
          typename ForwardIterator,
          typename PointPMap,
          typename NormalPMap>
OutputIterator
edge_aware_upsample_point_set(
  ForwardIterator first,    ///< forward iterator to the first input point.
  ForwardIterator beyond,   ///< past-the-end iterator.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap, ///< property map: `ForwardIterator` -> Point_3.
  NormalPMap normal_pmap, ///< property map: `ForwardIterator` -> Vector_3.
  double sharpness_angle,  ///< control sharpness(0-90)
  double edge_sensitivity,  ///< edge senstivity(0-5)
  double neighbor_radius, ///< initial size of neighbors.
  const unsigned int number_of_output_points///< number of iterations.   
  )
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return edge_aware_upsample_point_set(
    first, beyond,
    output,
    point_pmap,
    normal_pmap,
    sharpness_angle, 
    edge_sensitivity,
    neighbor_radius, 
    number_of_output_points,
    Kernel());
}
/// @endcond


/// @cond SKIP_IN_MANUAL
template <typename OutputIterator,
          typename ForwardIterator,
          typename NormalPMap>
OutputIterator
edge_aware_upsample_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  OutputIterator output, ///< output iterator over points.
  NormalPMap normal_pmap, ///< property map:  OutputIterator -> Vector_3.
  double sharpness_angle = 30,  ///< control sharpness(0-90)
  double edge_sensitivity = 1,  ///< edge senstivity(0-5)
  double neighbor_radius = -1, ///< initial size of neighbors.
  const unsigned int number_of_output_points = 1000///< number of output points.     
  )
{
  // just deduce value_type of OutputIterator
  return edge_aware_upsample_point_set
    <typename value_type_traits<OutputIterator>::type>(
    first, beyond,
    output,
    normal_pmap,
    sharpness_angle, 
    edge_sensitivity,
    neighbor_radius, 
    number_of_output_points);
}
/// @endcond

} //namespace CGAL

#endif // CGAL_UPSAMPLE_POINT_SET_H
