// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_PCA_ESTIMATE_NORMALS_H
#define CGAL_PCA_ESTIMATE_NORMALS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/IO/trace.h>
#include <CGAL/Dimension.h>
#include <CGAL/Point_set_processing_3/internal/Neighbor_query.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Memory_sizer.h>
#include <functional>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iterator>
#include <list>

#ifdef CGAL_LINKED_WITH_TBB
#include <CGAL/Point_set_processing_3/internal/Parallel_callback.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>  
#endif // CGAL_LINKED_WITH_TBB

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
template <typename NeighborQuery>
typename NeighborQuery::Kernel::Vector_3
pca_estimate_normal(const typename NeighborQuery::Kernel::Point_3& query, ///< point to compute the normal at
                    const NeighborQuery& neighbor_query, ///< KD-tree
                    unsigned int k, ///< number of neighbors
                    typename NeighborQuery::Kernel::FT neighbor_radius)
{
  // basic geometric types
  typedef typename NeighborQuery::Kernel Kernel;
  typedef typename Kernel::Point_3  Point;
  typedef typename Kernel::Plane_3  Plane;

  std::vector<Point> points;
  neighbor_query.get_points (query, k, neighbor_radius, std::back_inserter(points));
  
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

/**
   \ingroup PkgPointSetProcessing3Algorithms
   Estimates normal directions of the range of `points`
   by linear least squares fitting of a plane over the nearest neighbors.
   The output normals are randomly oriented.

   \pre `k >= 2`

   \tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are `Sequential_tag`,
                          `Parallel_tag`, and `Parallel_if_available_tag`.
   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param k number of neighbors
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin \cgalParamBegin{point_map} a model of
     `ReadablePropertyMap` with value type `geom_traits::Point_3`.  If
     this parameter is omitted,
     `CGAL::Identity_property_map<geom_traits::Point_3>` is
     used.\cgalParamEnd \cgalParamBegin{normal_map} a model of
     `WritablePropertyMap` with value type
     `geom_traits::Vector_3`.\cgalParamEnd
     \cgalParamBegin{neighbor_radius} spherical neighborhood
     radius. If provided, the neighborhood of a query point is
     computed with a fixed spherical radius instead of a fixed number
     of neighbors. In that case, the parameter `k` is used as a limit
     on the number of points returned by each spherical query (to
     avoid overly large number of points in high density areas). If no
     limit is wanted, use `k=0`.\cgalParamEnd
     \cgalParamBegin{callback} an instance of
     `std::function<bool(double)>`. It is called regularly when the
     algorithm is running: the current advancement (between 0. and 1.)
     is passed as parameter. If it returns `true`, then the algorithm
     continues its execution normally; if it returns `false`, the
     algorithm is stopped and the remaining normals are left
     unchanged.\cgalParamEnd \cgalParamBegin{geom_traits} an instance
     of a geometric traits class, model of `Kernel`\cgalParamEnd
     \cgalNamedParamsEnd
*/
template <typename ConcurrencyTag,
	  typename PointRange,
          typename NamedParameters
>
void
pca_estimate_normals(
  PointRange& points,
  unsigned int k,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_TRACE("Calls pca_estimate_normals()\n");

  // basic geometric types
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type NormalMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;
  typedef typename Kernel::FT FT;

  CGAL_static_assertion_msg(!(boost::is_same<NormalMap,
                              typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::NoMap>::value),
                            "Error: no normal map");

  PointMap point_map = choose_parameter(get_parameter(np, internal_np::point_map), PointMap());
  NormalMap normal_map = choose_parameter(get_parameter(np, internal_np::normal_map), NormalMap());
  FT neighbor_radius = choose_parameter(get_parameter(np, internal_np::neighbor_radius), FT(0));
  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                                 std::function<bool(double)>());

  typedef typename Kernel::Point_3 Point;

  // Input points types
  typedef typename PointRange::iterator iterator;
  typedef typename iterator::value_type value_type;
  typedef typename boost::property_traits<NormalMap>::value_type Vector;

  // types for K nearest neighbors search structure
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange, PointMap> Neighbor_query;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2);

  std::size_t memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("  Creates KD-tree\n");

  Neighbor_query neighbor_query (points, point_map);

  memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("  Computes normals\n");

  std::size_t nb_points = points.size();

  // iterate over input points, compute and output normal
  // vectors (already normalized)
#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
			     "Parallel_tag is enabled but TBB is unavailable.");
#else
  if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
  {
    Point_set_processing_3::internal::Parallel_callback
      parallel_callback (callback, nb_points);
     
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nb_points),
                      [&](const tbb::blocked_range<size_t>& r)
                      {
                        for( std::size_t i = r.begin(); i != r.end(); ++i)
                        {
                          if (parallel_callback.interrupted())
                            break;
                            
                          put (normal_map, *(points.begin() + i),
                               CGAL::internal::pca_estimate_normal
                               (get(point_map, *(points.begin() + i)), neighbor_query, k, neighbor_radius));
                            
                          ++ parallel_callback.advancement();
                        }
                      });

    parallel_callback.join();
  }
  else
#endif
  {
    std::size_t nb = 0;
    for (const value_type& vt : points)
    {
      put (normal_map, vt,
           internal::pca_estimate_normal(      
             get(point_map, vt),
             neighbor_query,
             k, neighbor_radius));
      if (callback && !callback ((nb+1) / double(nb_points)))
        break;
      ++ nb;
    }
  }
   
  memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("End of pca_estimate_normals()\n");
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename ConcurrencyTag,
	  typename PointRange
>
void
pca_estimate_normals(
  PointRange& points,
  unsigned int k) ///< number of neighbors.
{
  return pca_estimate_normals<ConcurrencyTag>
    (points, k, CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_PCA_ESTIMATE_NORMALS_H
