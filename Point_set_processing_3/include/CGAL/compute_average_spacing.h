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

#ifndef CGAL_AVERAGE_SPACING_3_H
#define CGAL_AVERAGE_SPACING_3_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Point_set_processing_3/internal/Neighbor_query.h>
#include <CGAL/Point_set_processing_3/internal/Callback_wrapper.h>
#include <CGAL/for_each.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/assertions.h>
#include <functional>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/iterator/zip_iterator.hpp>

#include <iterator>
#include <list>



#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#endif

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL
namespace internal {


/// Computes average spacing of one query point from K nearest neighbors.
///
/// \pre `k >= 2`.
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return average spacing (scalar).
template <typename NeighborQuery>
typename NeighborQuery::Kernel::FT
compute_average_spacing(const typename NeighborQuery::Kernel::Point_3& query, ///< 3D point whose spacing we want to compute
                        const NeighborQuery& neighbor_query,                      ///< KD-tree
                        unsigned int k)                        ///< number of neighbors
{
  // basic geometric types
  typedef typename NeighborQuery::Kernel Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;


  // performs k + 1 queries (if unique the query point is
  // output first). search may be aborted when k is greater
  // than number of input points
  FT sum_distances = (FT)0.0;
  unsigned int i = 0;
  neighbor_query.get_points
    (query, k, 0,
     boost::make_function_output_iterator
     ([&](const Point& p)
      {
        sum_distances += std::sqrt(CGAL::squared_distance (query,p));
        ++ i;
      }));

  // output average spacing
  return sum_distances / (FT)i;
}


} /* namespace internal */
/// \endcond



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessing3Algorithms
   Computes average spacing from k nearest neighbors.

   \pre `k >= 2.`

   \tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are `Sequential_tag`,
                          `Parallel_tag`, and `Parallel_if_available_tag`.
   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param k number of neighbors.
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
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
                       `false`, the algorithm is stopped, the average spacing value estimated
                       on the processed subset is returned.}
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

   \return average spacing (scalar). The return type `FT` is a number type. It is
   either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
   or the geometric traits class deduced from the point property map
   of `points`.
*/
template <typename ConcurrencyTag,
          typename PointRange,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS
>
#ifdef DOXYGEN_RUNNING
  FT
#else
  typename Point_set_processing_3::GetK<PointRange, CGAL_BGL_NP_CLASS>::Kernel::FT
#endif
compute_average_spacing(
  const PointRange& points,
  unsigned int k,
  const CGAL_BGL_NP_CLASS& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // basic geometric types
  typedef typename PointRange::const_iterator iterator;
  typedef typename CGAL::GetPointMap<PointRange, CGAL_BGL_NP_CLASS>::const_type PointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, CGAL_BGL_NP_CLASS>::Kernel Kernel;

  PointMap point_map = choose_parameter(get_parameter(np, internal_np::point_map), PointMap());
  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                                 std::function<bool(double)>());

  // types for K nearest neighbors search structure
  typedef typename Kernel::FT FT;
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, const PointRange&, PointMap> Neighbor_query;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2);

  // Instanciate a KD-tree search.
  Neighbor_query neighbor_query (points, point_map);

  // iterate over input points, compute and output normal
  // vectors (already normalized)
  FT sum_spacings = (FT)0.0;
  std::size_t nb = 0;
  std::size_t nb_points = std::distance(points.begin(), points.end());

  Point_set_processing_3::internal::Callback_wrapper<ConcurrencyTag>
    callback_wrapper (callback, nb_points);

  std::vector<FT> spacings (nb_points, -1);

  typedef boost::zip_iterator<boost::tuple<iterator, typename std::vector<FT>::iterator> > Zip_iterator;

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple (points.begin(), spacings.begin())),
                       boost::make_zip_iterator (boost::make_tuple (points.end(), spacings.end()))),
     [&](const typename Zip_iterator::reference& t)
     {
       if (callback_wrapper.interrupted())
         return false;

       get<1>(t) = CGAL::internal::compute_average_spacing<Neighbor_query>
         (get(point_map, get<0>(t)), neighbor_query, k);
       ++ callback_wrapper.advancement();

       return true;
     });

  for (unsigned int i = 0; i < spacings.size (); ++ i)
    if (spacings[i] >= 0.)
    {
      sum_spacings += spacings[i];
      ++ nb;
    }
  callback_wrapper.join();

  // return average spacing
  return sum_spacings / (FT)(nb);
}

/// \cond SKIP_IN_MANUAL

// variant with default NP
template <typename ConcurrencyTag, typename PointRange>
typename Point_set_processing_3::GetFT<PointRange>::type
compute_average_spacing(
  const PointRange& points,
  unsigned int k) ///< number of neighbors.
{
  return compute_average_spacing<ConcurrencyTag>
    (points, k, CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_AVERAGE_SPACING_3_H
