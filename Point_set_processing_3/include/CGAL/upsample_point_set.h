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
// Author(s) : Shihao Wu, Cl¨¦ment Jamin 

#ifndef CGAL_UPSAMPLE_POINT_SET_H
#define CGAL_UPSAMPLE_POINT_SET_H

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
namespace upsample_internal{


/// For each query point, select a best "base point" in its neighborhoods.
/// Then,a new point will be interpolated between query point and "base point".
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
  const rich_grid_internel::Rich_point<Kernel>& query, ///< 3D point to project
  const std::vector<rich_grid_internel::Rich_point<Kernel> >& 
                    neighbor_points,///< neighbor sample points
  const typename Kernel::FT edge_senstivity,
  unsigned int& output_base_index ///< base point index
  )
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;
  typedef typename rich_grid_internel::Rich_point<Kernel> Rich_point;


  FT best_dist2 = -10.0;
  Rich_point& v = query;
  for (unsigned int i = 0; i < neighbor_points.size(); i++)
  {
    Rich_point& t = neighbor_points[i];
    Point mid_point = (v.pt + t.pt) / FT(2.0);

    Vector& vm = v.normal;
    Vector& tm = t.normal;
    FT dot_produce = pow((FT(2.0) - vm * tm), edge_senstivity);

    Vector diff_t_mid = mid_point - t.pt;
    FT project_t = diff_t_mid * tm;
    FT min_dist2 = diff_t_mid.squared_length() - project_t * project_t;

    for (unsigned int j = 0; j < neighbor_points.size(); j++)
    {
      Rich_point& s = neighbor_points[j];
      Vector diff_s_mid = mid_point - s.pt;
      FT prject_s = diff_s_mid * s.normal;

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
      output_base_index = neighbor_points[i].index;
    }
  }

  return best_dist2;
}

} // namespace upsample_internal



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/// \ingroup PkgPointSetProcessing
/// Upsampling Algorithm:  progressively upsample the point set while 
/// approaching the edge singularities. 
/// More details please see: http://web.siat.ac.cn/~huihuang/EAR/EAR_page.html
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
template <typename ForwardIterator, 
          typename PointPMap, 
          typename NormalPMap,
          typename Kernel>
ForwardIterator
upsample_point_set(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  NormalPMap normal_pmap, ///< property map ForwardIterator -> Vector_3.
  const typename Kernel::FT sharpness_sigma,  ///< control sharpness(0-90)
  const typename Kernel::FT edge_senstivity,  ///< edge senstivity(0-5)
  const typename Kernel::FT neighbor_radius, ///< initial size of neighbors.
  const unsigned int number_of_output,///< number of iterations.                            
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
  typedef typename rich_grid_internel::Rich_box<Kernel> Rich_box;


  // preconditions
  CGAL_point_set_processing_precondition(first != beyond);
  CGAL_point_set_processing_precondition(sharpness_sigma >= 0 
                                       &&sharpness_sigma <= 90);
  CGAL_point_set_processing_precondition(edge_senstivity >= 0 
                                       &&edge_senstivity <= 5);
  CGAL_point_set_processing_precondition(neighbor_radius > 0);

  std::size_t number_of_input = std::distance(first, beyond);
  CGAL_point_set_processing_precondition(number_of_output > number_of_input);

  // copy rich point set
  ForwardIterator it;// point iterator
  unsigned int i; 
  std::vector<Rich_point> rich_point_set(number_of_input);
  Rich_box box;
  for(it = first, i = 0; it != beyond; ++it, i++)
  {
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    rich_point_set[i].pt = get(point_pmap, it);
    rich_point_set[i].normal = get(normal_pmap, it);
#else
    rich_point_set[i].pt = get(point_pmap, *it);
    rich_point_set[i].normal = get(normal_pmap, *it);
#endif

    rich_point_set[i].index = i;
    box.add_point(rich_point_set[i].pt);
  }

  // compute neighborhood
  rich_grid_internel::compute_ball_neighbors_one_self(rich_point_set,
                                                      box,
                                                      neighbor_radius);





  return first;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename ForwardIterator, typename PointPMap, typename NormalPMap>
ForwardIterator
upsample_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  NormalPMap normal_pmap, ///< property map ForwardIterator -> Vector_3.
  double sharpness_sigma,  ///< control sharpness(0-90)
  double edge_senstivity,  ///< edge senstivity(0-5)
  double neighbor_radius, ///< initial size of neighbors.
  const unsigned int number_of_output_points///< number of iterations.                                  
) 
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return upsample_point_set(
    first, beyond,
    point_pmap,
    normal_pmap,
    sharpness_sigma,
    edge_senstivity,
    neighbor_radius,
    number_of_output_points,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Dereference_property_map.
template <typename ForwardIterator, typename NormalPMap>
ForwardIterator
upsample_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  NormalPMap normal_pmap, ///< property map ForwardIterator -> Vector_3.
  double sharpness_sigma,  ///< control sharpness(0-90)
  double edge_senstivity,  ///< edge senstivity(0-5)
  double neighbor_radius, ///< initial size of neighbors.
  const unsigned int number_of_output_points///< number of iterations.   
)
{
  return upsample_point_set(
    first, beyond,
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    make_dereference_property_map(first),
#else
    make_identity_property_map(typename std::iterator_traits<ForwardIterator>::
                               value_type()),
#endif
    normal_pmap,
    sharpness_sigma, 
    edge_senstivity,
    neighbor_radius, 
    number_of_output_points
    );
}
/// @endcond

} //namespace CGAL

#endif // CGAL_UPSAMPLE_POINT_SET_H
