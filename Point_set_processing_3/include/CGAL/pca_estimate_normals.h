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
#include <CGAL/Point_set_processing_3/internal/Callback_wrapper.h>
#include <CGAL/for_each.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Memory_sizer.h>
#include <functional>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iterator>
#include <list>

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

   \param points input point range
   \param k number of neighbors
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point set `points`}
       \cgalParamType{a model of `WritablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Vector_3`}
     \cgalParamNEnd

     \cgalParamNBegin{neighbor_radius}
       \cgalParamDescription{the spherical neighborhood radius}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{`0` (no limit)}
       \cgalParamExtra{If provided, the neighborhood of a query point is computed with a fixed spherical
                       radius instead of a fixed number of neighbors. In that case, the parameter
                       `k` is used as a limit on the number of points returned by each spherical
                       query (to avoid overly large number of points in high density areas).}
     \cgalParamNEnd

     \cgalParamNBegin{callback}
       \cgalParamDescription{a mechanism to get feedback on the advancement of the algorithm
                             while it's running and to interrupt it if needed}
       \cgalParamType{an instance of `std::function<bool(double)>`.}
       \cgalParamDefault{unused}
       \cgalParamExtra{It is called regularly when the
                       algorithm is running: the current advancement (between 0. and
                       1.) is passed as parameter. If it returns `true`, then the
                       algorithm continues its execution normally; if it returns
                       `false`, the algorithm is stopped and the remaining normals are left unchanged.}
       \cgalParamExtra{The callback will be copied and therefore needs to be lightweight.}
       \cgalParamExtra{When `CGAL::Parallel_tag` is used, the `callback` mechanism is called asynchronously
                       on a separate thread and shouldn't access or modify the variables that are parameters of the algorithm.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd
*/
template <typename ConcurrencyTag,
          typename PointRange,
          typename NamedParameters = parameters::Default_named_parameters
>
void
pca_estimate_normals(
  PointRange& points,
  unsigned int k,
  const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_TRACE_STREAM << "Calls pca_estimate_normals()\n";

  // basic geometric types
  typedef Point_set_processing_3_np_helper<PointRange, NamedParameters> NP_helper;
  typedef typename NP_helper::Point_map PointMap;
  typedef typename NP_helper::Normal_map NormalMap;
  typedef typename NP_helper::Geom_traits Kernel;
  typedef typename Kernel::FT FT;

  CGAL_static_assertion_msg(NP_helper::has_normal_map(), "Error: no normal map");

  PointMap point_map = NP_helper::get_point_map(points, np);
  NormalMap normal_map = NP_helper::get_normal_map(points, np);
  FT neighbor_radius = choose_parameter(get_parameter(np, internal_np::neighbor_radius), FT(0));
  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                                 std::function<bool(double)>());

  // Input points types
  typedef typename PointRange::iterator iterator;
  typedef typename iterator::value_type value_type;

  // types for K nearest neighbors search structure
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange&, PointMap> Neighbor_query;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2);

  std::size_t memory = CGAL::Memory_sizer().virtual_size();
  CGAL_TRACE_STREAM << (memory >> 20) << " Mb allocated\n";
  CGAL_TRACE_STREAM << "  Creates KD-tree\n";

  Neighbor_query neighbor_query (points, point_map);

  memory = CGAL::Memory_sizer().virtual_size();
  CGAL_TRACE_STREAM << (memory >> 20) << " Mb allocated\n";
  CGAL_TRACE_STREAM << "  Computes normals\n";

  std::size_t nb_points = points.size();

  Point_set_processing_3::internal::Callback_wrapper<ConcurrencyTag>
    callback_wrapper (callback, nb_points);

  CGAL::for_each<ConcurrencyTag>
    (points,
     [&](value_type& vt)
     {
       if (callback_wrapper.interrupted())
         return false;

       put (normal_map, vt,
            CGAL::internal::pca_estimate_normal
            (get(point_map, vt), neighbor_query, k, neighbor_radius));

       ++ callback_wrapper.advancement();

       return true;
     });

  callback_wrapper.join();

  memory = CGAL::Memory_sizer().virtual_size();
  CGAL_TRACE_STREAM << (memory >> 20) << " Mb allocated\n";
  CGAL_TRACE_STREAM << "End of pca_estimate_normals()\n";
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_PCA_ESTIMATE_NORMALS_H
