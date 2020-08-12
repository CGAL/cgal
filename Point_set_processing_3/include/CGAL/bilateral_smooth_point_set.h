// Copyright (c) 2013-06  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Shihao Wu, Clement Jamin, Pierre Alliez

#ifndef CGAL_BILATERAL_SMOOTH_POINT_SET_H
#define CGAL_BILATERAL_SMOOTH_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/number_type_config.h>
#include <CGAL/Point_set_processing_3/internal/Neighbor_query.h>
#include <CGAL/Point_set_processing_3/internal/Callback_wrapper.h>
#include <CGAL/for_each.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/squared_distance_3.h>
#include <functional>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/iterator/zip_iterator.hpp>

#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <CGAL/Real_timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/property_map.h>

//#define CGAL_PSP3_VERBOSE

namespace CGAL {

// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL

namespace bilateral_smooth_point_set_internal{

/// Compute bilateral projection for each point
/// according to their KNN neighborhood points
///
/// \pre `k >= 2`, radius > 0 , sharpness_angle > 0 && sharpness_angle < 90
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return

template <typename Kernel, typename PointRange,
          typename PointMap, typename VectorMap>
std::pair<typename Kernel::Point_3, typename Kernel::Vector_3>
compute_denoise_projection(
  const typename PointRange::iterator::value_type& vt,
  PointMap point_map,
  VectorMap normal_map,
  const std::vector<typename PointRange::iterator>& neighbor_pwns,
  typename Kernel::FT radius,                   ///< accept neighborhood radius
  typename Kernel::FT sharpness_angle           ///< control sharpness(0-90)
)
{
  CGAL_point_set_processing_precondition(radius > 0);
  CGAL_point_set_processing_precondition(sharpness_angle > 0
                                         && sharpness_angle < 90);

  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Point_3 Point;

  FT radius2 = radius * radius;

  FT weight = (FT)0.0;
  FT iradius16 = -(FT)4.0/radius2;
  FT project_dist_sum = FT(0.0);
  FT project_weight_sum = FT(0.0);
  Vector normal_sum = CGAL::NULL_VECTOR;

  FT cos_sigma = cos(sharpness_angle * CGAL_PI / 180.0);
  FT sharpness_bandwidth = std::pow((CGAL::max)(1e-8, 1 - cos_sigma), 2);

  for (typename PointRange::iterator it : neighbor_pwns)
  {
    const Point& np = get(point_map, *it);
    const Vector& nn = get(normal_map, *it);

    FT dist2 = CGAL::squared_distance(get(point_map, vt), np);
    if (dist2 < radius2)
    {
      FT theta = std::exp(dist2 * iradius16);
      FT psi = std::exp(-std::pow(1 - get(normal_map, vt) * nn, 2)
        / sharpness_bandwidth);

      weight = theta * psi;

      project_dist_sum += ((get(point_map, vt) - np) * nn) * weight;
      project_weight_sum += weight;
      normal_sum = normal_sum + nn * weight;
    }
  }

  Vector update_normal = normal_sum / project_weight_sum;
  update_normal = update_normal / sqrt(update_normal.squared_length());

  Point update_point = get(point_map, vt) - update_normal *
                      (project_dist_sum / project_weight_sum);

  return std::make_pair (update_point, update_normal);
}

/// Computes max-spacing of one query point from K nearest neighbors.
///
/// \pre `k >= 2`.
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return max spacing.
template <typename NeighborQuery>
typename NeighborQuery::Kernel::FT
compute_max_spacing(
  const typename NeighborQuery::value_type& vt,
  typename NeighborQuery::Point_map point_map,
  NeighborQuery& neighbor_query,                                     ///< KD-tree
  unsigned int k)                                 ///< number of neighbors
{
  // basic geometric types
  typedef typename NeighborQuery::Kernel Kernel;
  typedef typename Kernel::FT FT;

  // performs k + 1 queries (if unique the query point is
  // output first). search may be aborted when k is greater
  // than number of input points
  FT max_distance = (FT)0.0;
  neighbor_query.get_iterators
    (get(point_map, vt), k, (FT)(0.0),
     boost::make_function_output_iterator
     ([&](const typename NeighborQuery::input_iterator& it)
      {
        double dist2 = CGAL::squared_distance (get(point_map, vt), get(point_map, *it));
        max_distance = (CGAL::max)(dist2, max_distance);
      }));

  // output max spacing
  return std::sqrt(max_distance);
}

} /* namespace internal */

/// \endcond


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/**
   \ingroup PkgPointSetProcessing3Algorithms

   This function smooths an input point set by iteratively projecting each
   point onto the implicit surface patch fitted over its nearest neighbors.
   Bilateral projection preserves sharp features according to the normal
   (gradient) information. Both point positions and normals will be modified.
   For more details, please see section 4 in \cgalCite{ear-2013}.

   A parallel version of this function is provided and requires the executable to be
   linked against the <a href="https://www.threadingbuildingblocks.org">Intel TBB library</a>.
   To control the number of threads used, the user may use the tbb::task_scheduler_init class.
   See the <a href="https://www.threadingbuildingblocks.org/documentation">TBB documentation</a>
   for more details.

   \pre Normals must be unit vectors
   \pre k >= 2

   \tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are `Sequential_tag`,
                          `Parallel_tag`, and `Parallel_if_available_tag`.
   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param k size of the neighborhood for the implicit surface patch fitting.
   The larger the value is, the smoother the result will be.
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadWritePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point set `points`}
       \cgalParamType{a model of `ReadWritePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Vector_3`}
       \cgalParamDefault{Normals are computed and stored internally.}
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

     \cgalParamNBegin{sharpness_angle}
       \cgalParamDescription{controls the sharpness of the result}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{`30`}
       \cgalParamExtra{The larger the value is, the smoother the result will be.
                       The range of possible value is `[0, 90]`}
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
                       `false`, the algorithm is stopped, all points are left unchanged
                       and the function return `NaN`.}
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

   \return Average point movement error. It's a convergence criterium for the algorithm.
   This value can help the user to decide how many iterations are
   sufficient.
*/
template <typename ConcurrencyTag,
          typename PointRange,
          typename NamedParameters>
double
bilateral_smooth_point_set(
  PointRange& points,
  unsigned int k,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // basic geometric types
  typedef typename PointRange::iterator iterator;
  typedef typename iterator::value_type value_type;
  typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type NormalMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;

  CGAL_static_assertion_msg(!(boost::is_same<NormalMap,
                              typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::NoMap>::value),
                            "Error: no normal map");

  typedef typename Kernel::FT FT;

  double sharpness_angle = choose_parameter(get_parameter(np, internal_np::sharpness_angle), 30.);
  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                                 std::function<bool(double)>());

  CGAL_point_set_processing_precondition(points.begin() != points.end());
  CGAL_point_set_processing_precondition(k > 1);

  // types for K nearest neighbors search structure
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange&, PointMap> Neighbor_query;

  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  NormalMap normal_map = choose_parameter<NormalMap>(get_parameter(np, internal_np::normal_map));
  FT neighbor_radius = choose_parameter(get_parameter(np, internal_np::neighbor_radius), FT(0));

  std::size_t nb_points = points.size();

#ifdef CGAL_PSP3_VERBOSE
   std::cout << "Initialization and compute max spacing: " << std::endl;
#endif
   // initiate a KD-tree search for points
   Neighbor_query neighbor_query (points, point_map);

   // Guess spacing
#ifdef CGAL_PSP3_VERBOSE
   CGAL::Real_timer task_timer;
   task_timer.start();
#endif
   FT guess_neighbor_radius = 0.0;

   for (const value_type& vt : points)
   {
     FT max_spacing = bilateral_smooth_point_set_internal::
       compute_max_spacing (vt, point_map, neighbor_query, k);
     guess_neighbor_radius = (CGAL::max)(max_spacing, guess_neighbor_radius);
   }

#ifdef CGAL_PSP3_VERBOSE
   task_timer.stop();
#endif
   guess_neighbor_radius *= 0.95;

#ifdef CGAL_PSP3_VERBOSE
   CGAL::Memory_sizer::size_type memory = CGAL::Memory_sizer().virtual_size();
   std::cout << "done: " << task_timer.time() << " seconds, "
             << (memory>>20) << " Mb allocated" << std::endl;

   std::cout << "Compute all neighbors: " << std::endl;
   task_timer.reset();
   task_timer.start();
#endif
   // compute all neighbors
   typedef std::vector<iterator> iterators;
   std::vector<iterators> pwns_neighbors;
   pwns_neighbors.resize(nb_points);

   Point_set_processing_3::internal::Callback_wrapper<ConcurrencyTag>
     callback_wrapper (callback, 2 * nb_points);

   typedef boost::zip_iterator<boost::tuple<iterator, typename std::vector<iterators>::iterator> > Zip_iterator;

   CGAL::for_each<ConcurrencyTag>
     (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple (points.begin(), pwns_neighbors.begin())),
                        boost::make_zip_iterator (boost::make_tuple (points.end(), pwns_neighbors.end()))),
      [&](const typename Zip_iterator::reference& t)
      {
        if (callback_wrapper.interrupted())
          return false;

        neighbor_query.get_iterators (get(point_map, get<0>(t)), k, neighbor_radius,
                                      std::back_inserter (get<1>(t)));

        ++ callback_wrapper.advancement();

        return true;
      });

   bool interrupted = callback_wrapper.interrupted();

   // We interrupt by hand as counter only goes halfway and won't terminate by itself
   callback_wrapper.interrupted() = true;
   callback_wrapper.join();

   // If interrupted during this step, nothing is computed, we return NaN
   if (interrupted)
     return std::numeric_limits<double>::quiet_NaN();

#ifdef CGAL_PSP3_VERBOSE
   task_timer.stop();
   memory = CGAL::Memory_sizer().virtual_size();
   std::cout << "done: " << task_timer.time() << " seconds, "
             << (memory>>20) << " Mb allocated" << std::endl;

   std::cout << "Compute update points and normals: " << std::endl;
   task_timer.reset();
   task_timer.start();
#endif
   // update points and normals
   std::vector<std::pair<Point_3, Vector_3> > update_pwns(nb_points);

   callback_wrapper.reset (2 * nb_points, nb_points);

   typedef boost::zip_iterator
     <boost::tuple<iterator,
                   typename std::vector<iterators>::iterator,
                   typename std::vector<std::pair<Point_3, Vector_3> >::iterator> > Zip_iterator_2;


   CGAL::for_each<ConcurrencyTag>
     (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple
                                                  (points.begin(), pwns_neighbors.begin(), update_pwns.begin())),
                        boost::make_zip_iterator (boost::make_tuple
                                                  (points.end(), pwns_neighbors.end(), update_pwns.end()))),
      [&](const typename Zip_iterator_2::reference& t)
      {
        if (callback_wrapper.interrupted())
          return false;

        get<2>(t) = bilateral_smooth_point_set_internal::
          compute_denoise_projection<Kernel, PointRange>
          (get<0>(t),
           point_map, normal_map,
           get<1>(t),
           guess_neighbor_radius,
           sharpness_angle);

        ++ callback_wrapper.advancement();

        return true;
      });

   callback_wrapper.join();

   // If interrupted during this step, nothing is computed, we return NaN
   if (callback_wrapper.interrupted())
     return std::numeric_limits<double>::quiet_NaN();

#ifdef CGAL_PSP3_VERBOSE
   task_timer.stop();
   memory = CGAL::Memory_sizer().virtual_size();
   std::cout << "done: " << task_timer.time() << " seconds, "
             << (memory>>20) << " Mb allocated" << std::endl;
#endif
   // save results
   FT sum_move_error = 0;
   std::size_t nb = 0;
   for (value_type& vt : points)
   {
     sum_move_error += CGAL::squared_distance(get(point_map, vt), update_pwns[nb].first);
     put (point_map, vt, update_pwns[nb].first);
     put (normal_map, vt, update_pwns[nb].second);
     ++ nb;
   }

   return sum_move_error / nb_points;
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename ConcurrencyTag,
          typename PointRange>
double
bilateral_smooth_point_set(
  PointRange& points,
  unsigned int k)           ///< size of the neighborhood for the implicit surface patch fitting.
                            ///< The larger the value is, the smoother the result will be.
{
  return bilateral_smooth_point_set<ConcurrencyTag>
    (points, k, CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_BILATERAL_SMOOTH_POINT_SET_H
