// Copyright (c) 2019  GeometryFactory(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Necip Fazil Yildiran

#ifndef CGAL_POINTMATCHER_REGISTER_POINT_SETS_H
#define CGAL_POINTMATCHER_REGISTER_POINT_SETS_H

#include <CGAL/license/Point_set_processing_3.h>

#ifdef CGAL_LINKED_WITH_POINTMATCHER

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/pointmatcher/compute_registration_transformation.h>

#include <boost/type_traits/is_same.hpp>

#include <Eigen/Dense>

namespace CGAL {

namespace pointmatcher {

// point_set_1 is reference while point_set_2 is data
/**
   \ingroup PkgPointSetProcessing3Algorithms
   
   Computes the registration of `point_set_2` with respect to `point_set_1` and
   applies it.

   Registration is computed using the Iterative Closest Point (ICP) algorithm.

   \note This function requires the \ref thirdpartylibpointmatcher library.

   \tparam PointRange1 is a model of `Range`. The value type of its iterator is
   the key type of the named parameter `point_map` in `NamedParameters1`.
   \tparam PointRange2 is a model of `Range`. The value type of its iterator is
   the key type of the named parameter `point_map` in `NamedParameters2`.

   \param point_set_1 input point range used as reference.
   \param point_set_2 input point range whose registration w.r.t. `point_set_1` will be computed.
   \param np1 optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` whose key type
     is the value type of the iterator of `PointRange1` and whose value type is
     `geom_traits::Point_3`.  If this parameter is omitted,
     `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{normal_map} a model of `ReadablePropertyMap` whose key
     type is the value type of the iterator of `PointRange1` and whose value
     type `geom_traits::Vector_3`.\cgalParamEnd

     \cgalParamBegin{point_set_filters} is a model of `Range`. The value type of 
     its iterator is `ICP_config`. 
     
     The chain of filters to be applied to the reference point cloud. The reference
     point cloud is processed into an intermediate point cloud with the given chain
     of filters to be used in the alignment procedure. The chain is organized with
     the forward traversal order of the point set filters range.
     
     The chain of point set filters are applied only once at the beginning of the
     ICP procedure, i.e., before the first iteration of the ICP algorithm.
     
     The filters can have several purposes, including but are not limited to
     i) removal of noisy points which render alignment of point clouds difficult,
     ii) removal of redundant points so as to speed up alignment, iii) addition 
     of descriptive information to the points such as a surface normal vector, 
     or the direction from the point to the sensor.

     Corresponds to `referenceDataPointsFilters` configuration module of \ref thirdpartylibpointmatcher
     library. The filters should be chosen and set from possible components of
     the `referenceDataPointsFilters` configuration module. 
     See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
     for possible configurations.

     If this parameter is omitted, `RandomSamplingDataPointsFilter` is used.
     \cgalParamEnd

     \cgalParamBegin{matcher} is a model of `ICP_config`.
     The method used for matching (linking) the points from `point_set_2`, to the 
     points in the reference cloud, `point_set_1`.
     
     Corresponds to `matcher` configuration module of \ref thirdpartylibpointmatcher
     library. The matcher should be chosen and set from possible components of
     the `matcher` configuration module. 
     See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
     for possible configurations.

     If this parameter is omitted, `KDTreeMatcher` is used.
     \cgalParamEnd

     \cgalParamBegin{outlier_filters} is a model of `Range`. The value type of 
     its iterator is `ICP_config`.
     The chain of filters to be applied to the matched (linked) point clouds after
     each processing iteration of the ICP algorithm to remove the links which do not
     correspond to true point correspondences. The outliers are rejected. Points
     with no link are ignored in the subsequent error minimization step.
     The chain is organized with the forward traversal order of the outlier filters
     range.

     Corresponds to `outlierFilters` configuration module of \ref thirdpartylibpointmatcher
     library. The filters should be chosen and set from possible components of
     the `outlierFilters` configuration module. 
     See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
     for possible configurations.

     If this parameter is omitted, `TrimmedDistOutlierFilter` is used.
     \cgalParamEnd

     \cgalParamBegin{error_minimizer} is a model of `ICP_config`.
     The error minimizer that computes a transformation matrix such as to minimize
     the error between the point sets.

     Corresponds to `errorMinimizer` configuration module of \ref thirdpartylibpointmatcher
     library. The error minimizer should be chosen and set from possible components of
     the `errorMinimizer` configuration module. 
     See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
     for possible configurations.

     If this parameter is omitted, `PointToPlaneErrorMinimizer` is used.
     \cgalParamEnd

     \cgalParamBegin{transformation_checkers} is a model of `Range`. The value type of 
     its iterator is `ICP_config`.
     The chain of transformation checkers. A transformation checker can stop the
     iteration depending on the conditions it defines.
     
     The chain is organized with the forward traversal order of the transformation
     checkers range.

     Corresponds to `transformationCheckers` configuration module of \ref thirdpartylibpointmatcher
     library. The transformation checkers should be chosen and set from possible components of
     the `transformationCheckers` configuration module. 
     See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
     for possible configurations.

     If this parameter is omitted, the chain of `CounterTransformationChecker` and
     `DifferentialTransformationChecker` is used.
     \cgalParamEnd
     
     \cgalParamBegin{inspector} is a model of `ICP_config`.
     The inspector allows to log data at different steps for analysis. Inspectors
     typically provide deeper scrutiny than the logger.

     Corresponds to `inspector` configuration module of \ref thirdpartylibpointmatcher
     library. The inspector should be chosen and set from possible components of
     the `inspector` configuration module. 
     See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
     for possible configurations.

     If this parameter is omitted, `NullInspector` is used.
     \cgalParamEnd

     \cgalParamBegin{logger} is a model of `ICP_config`.
     The method for logging information regarding the registration process outputted
     by \ref thirdpartylibpointmatcher library. The logs generated by CGAL library
     does not get effected by this configuration.

     Corresponds to `logger` configuration module of \ref thirdpartylibpointmatcher
     library. The logger should be chosen and set from possible components of
     the `logger` configuration module. 
     See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
     for possible configurations.

     If this parameter is omitted, `NullLogger` is used.
     \cgalParamEnd

     \cgalParamBegin{geom_traits} an instance of a geometric traits class,
     model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

   \param np2 optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.
   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` whose key type
     is the value type of the iterator of `PointRange2` and whose value type is
     `geom_traits::Point_3`.  If this parameter is omitted,
     `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{normal_map} a model of `ReadablePropertyMap` whose key
     type is the value type of the iterator of `PointRange2` and whose value
     type `geom_traits::Vector_3`.\cgalParamEnd

     \cgalParamBegin{point_set_filters} is a model of `Range`. The value type of 
     its iterator is `ICP_config`.

     The chain of filters to be applied to the point cloud, `point_set_2`. The
     point cloud is processed into an intermediate point cloud with the given chain
     of filters to be used in the alignment procedure. The chain is organized with
     the forward traversal order of the point set filters range.
     
     The chain of point set filters are applied only once at the beginning of the
     ICP procedure, i.e., before the first iteration of ICP algorithm.
     
     The filters can have several purposes, including but are not limited to
     i) removal of noisy points which render alignment of point clouds difficult,
     ii) removal of redundant points so as to speed up alignment, iii) addition 
     of descriptive information to the points such as a surface normal vector, 
     or the direction from the point to the sensor.

     Corresponds to `readingDataPointsFilters` configuration module of \ref thirdpartylibpointmatcher
     library. The filters should be chosen and set from possible components of
     the `readingDataPointsFilters` configuration module. 
     See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
     for possible configurations.

     If this parameter is omitted, `SamplingSurfaceNormalDataPointsFilter` is used.
     \cgalParamEnd

     \cgalParamBegin{transformation} The affine transformation that is used as the 
     initial transformation for `point_set_2`.

     If this parameter is omitted, identity transformation is used.
     \cgalParamEnd

   \cgalNamedParamsEnd
   \return `true` if registration is converged, `false` otherwise. A log why it
   failed to converge is written to `std::cerr` if the registration cannot converge. 
*/
template <class PointRange1, class PointRange2,
          class NamedParameters1, class NamedParameters2>
bool
register_point_sets (const PointRange1& point_set_1, PointRange2& point_set_2,
                     const NamedParameters1& np1, const NamedParameters2& np2)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  namespace PSP = CGAL::Point_set_processing_3;
  typedef typename PSP::GetK<PointRange1, NamedParameters1>::Kernel Kernel;

  // compute registration transformation
  std::pair<typename Kernel::Aff_transformation_3, bool> res =
    compute_registration_transformation(point_set_1, point_set_2, np1, np2);

  // property map type of point_set_2 
  typedef typename CGAL::GetPointMap<PointRange2, NamedParameters2>::type PointMap2;
  PointMap2 point_map2 = choose_parameter(get_parameter(np2, internal_np::point_map), PointMap2());

  // update CGAL points
  for (typename PointRange2::iterator it=point_set_2.begin(),
                                      end=point_set_2.end(); it!=end; ++it)
  {
    put(point_map2, *it, get(point_map2, *it).transform(res.first));
  }

  return res.second;
}

// convenience overloads
template <class PointRange1, class PointRange2,
          class NamedParameters1>
bool
register_point_sets(const PointRange1& point_set_1, PointRange2& point_set_2,
                    const NamedParameters1& np1)
{
  namespace params = CGAL::Point_set_processing_3::parameters;
  return register_point_sets(point_set_1, point_set_2, np1, params::all_default(point_set_1));
}

template <class PointRange1, class PointRange2>
bool
register_point_sets(const PointRange1& point_set_1, PointRange2& point_set_2)
{
  namespace params = CGAL::Point_set_processing_3::parameters;
  return register_point_sets(point_set_1, point_set_2,
                             params::all_default(point_set_1),
                             params::all_default(point_set_2));
}

} } // end of namespace CGAL::pointmatcher

#endif // CGAL_LINKED_WITH_POINTMATCHER

#endif // CGAL_POINTMATCHER_REGISTER_POINT_SETS_H
