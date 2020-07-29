// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez, Marc Pouget and Laurent Saboret

#ifndef CGAL_JET_SMOOTH_POINT_SET_H
#define CGAL_JET_SMOOTH_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/IO/trace.h>
#include <CGAL/Point_set_processing_3/internal/Neighbor_query.h>
#include <CGAL/Point_set_processing_3/internal/Callback_wrapper.h>
#include <CGAL/for_each.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <functional>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/iterator/zip_iterator.hpp>

#include <iterator>
#include <list>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL

namespace internal {

/// Smoothes one point position using jet fitting on the k
/// nearest neighbors and reprojection onto the jet.
///
/// \pre `k >= 2`
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return computed point
template <typename SvdTraits,
          typename NeighborQuery
          >
typename NeighborQuery::Kernel::Point_3
jet_smooth_point(
  const typename NeighborQuery::Kernel::Point_3& query, ///< 3D point to project
  NeighborQuery& neighbor_query, ///< KD-tree
  const unsigned int k, ///< number of neighbors.
  typename NeighborQuery::Kernel::FT neighbor_radius,
  const unsigned int degree_fitting,
  const unsigned int degree_monge)
{
  // basic geometric types
  typedef typename NeighborQuery::Kernel Kernel;
  typedef typename Kernel::Point_3 Point;

  // types for jet fitting
  typedef Monge_via_jet_fitting< Kernel,
                                 Simple_cartesian<double>,
                                 SvdTraits> Monge_jet_fitting;
  typedef typename Monge_jet_fitting::Monge_form Monge_form;

  std::vector<Point> points;

  // query using as fallback minimum requires nb points for jet fitting (d+1)*(d+2)/2
  neighbor_query.get_points (query, k, neighbor_radius, std::back_inserter(points),
                             (degree_fitting + 1) * (degree_fitting + 2) / 2);

  // performs jet fitting
  Monge_jet_fitting monge_fit;
  Monge_form monge_form = monge_fit(points.begin(), points.end(),
                                    degree_fitting, degree_monge);

  // output projection of query point onto the jet
  return monge_form.origin();
}


} /* namespace internal */

/// \endcond



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessing3Algorithms
   Smoothes the range of `points` using jet fitting on the
   nearest neighbors and reprojection onto the jet.
   As this method relocates the points, it
   should not be called on containers sorted w.r.t. point locations.

   \pre `k >= 2`

   \tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are `Sequential_tag`,
                          `Parallel_tag`, and `Parallel_if_available_tag`.
   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param k number of neighbors
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
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

     \cgalParamNBegin{degree_fitting}
       \cgalParamDescription{the degree of fitting}
       \cgalParamType{unsigned int}
       \cgalParamDefault{`2`}
       \cgalParamExtra{see `CGAL::Monge_via_jet_fitting`}
     \cgalParamNEnd

     \cgalParamNBegin{degree_monge}
       \cgalParamDescription{the Monge degree}
       \cgalParamType{unsigned int}
       \cgalParamDefault{`2`}
       \cgalParamExtra{see `CGAL::Monge_via_jet_fitting`}
     \cgalParamNEnd

     \cgalParamNBegin{svd_traits}
       \cgalParamDescription{the linear algebra algorithm used in the class `CGAL::Monge_via_jet_fitting`}
       \cgalParamType{a class fitting the requirements of `CGAL::Monge_via_jet_fitting`}
       \cgalParamDefault{If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available
                         and `CGAL_EIGEN3_ENABLED` is defined, then `CGAL::Eigen_svd` is used.}
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
                       `false`, the algorithm is stopped and the remaining points are left unchanged.}
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
          typename NamedParameters
>
void
jet_smooth_point_set(
  PointRange& points,
  unsigned int k,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // basic geometric types
  typedef typename PointRange::iterator iterator;
  typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;
  typedef typename GetSvdTraits<NamedParameters>::type SvdTraits;

  CGAL_static_assertion_msg(!(boost::is_same<SvdTraits,
                              typename GetSvdTraits<NamedParameters>::NoTraits>::value),
                            "Error: no SVD traits");

  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  typename Kernel::FT neighbor_radius = choose_parameter(get_parameter(np, internal_np::neighbor_radius),
                                                         typename Kernel::FT(0));
  unsigned int degree_fitting = choose_parameter(get_parameter(np, internal_np::degree_fitting), 2);
  unsigned int degree_monge = choose_parameter(get_parameter(np, internal_np::degree_monge), 2);
  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                               std::function<bool(double)>());

  // types for K nearest neighbors search structure
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange&, PointMap> Neighbor_query;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2);

  // Instanciate a KD-tree search.
  Neighbor_query neighbor_query (points, point_map);

  // Iterates over input points and mutates them.
  // Implementation note: the cast to Point& allows to modify only the point's position.

  std::size_t nb_points = points.size();

  Point_set_processing_3::internal::Callback_wrapper<ConcurrencyTag>
    callback_wrapper (callback, nb_points);

  std::vector<typename Kernel::Point_3> smoothed (points.size());

  typedef boost::zip_iterator
     <boost::tuple<iterator,
                   typename std::vector<typename Kernel::Point_3>::iterator> > Zip_iterator;

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple (points.begin(), smoothed.begin())),
                       boost::make_zip_iterator (boost::make_tuple (points.end(), smoothed.end()))),
     [&](const typename Zip_iterator::reference& t)
     {
       if (callback_wrapper.interrupted())
         return false;

       get<1>(t) = CGAL::internal::jet_smooth_point<SvdTraits>
         (get (point_map, get<0>(t)), neighbor_query,
          k,
          neighbor_radius,
          degree_fitting,
          degree_monge);
       ++ callback_wrapper.advancement();

       return true;
     });

  callback_wrapper.join();

  // Finally, update points
  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple (points.begin(), smoothed.begin())),
                       boost::make_zip_iterator (boost::make_tuple (points.end(), smoothed.end()))),
     [&](const typename Zip_iterator::reference& t)
     {
       put (point_map, get<0>(t), get<1>(t));
       return true;
     });
}


/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename ConcurrencyTag,
          typename PointRange>
void
jet_smooth_point_set(
  PointRange& points,
  unsigned int k) ///< number of neighbors.
{
  jet_smooth_point_set<ConcurrencyTag>
    (points, k, CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_JET_SMOOTH_POINT_SET_H
