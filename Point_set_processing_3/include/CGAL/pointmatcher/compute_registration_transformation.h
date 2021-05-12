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

#ifndef CGAL_POINTMATCHER_COMPUTE_REGISTRATION_TRANSFORMATION_H
#define CGAL_POINTMATCHER_COMPUTE_REGISTRATION_TRANSFORMATION_H

#include <CGAL/license/Point_set_processing_3.h>

#if defined(CGAL_LINKED_WITH_POINTMATCHER) || defined(DOXYGEN_RUNNING)

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/aff_transformation_tags.h>

#include <boost/type_traits/is_same.hpp>

#include <pointmatcher/PointMatcher.h>

#include <iostream>
#include <string>
#include <map>

namespace CGAL {

namespace pointmatcher {

template<typename Scalar>
using ICP = typename PointMatcher<Scalar>::ICP;

/*!
   \ingroup PkgPointSetProcessing3Algorithms

   \brief The class `ICP_config` is designed to handle preparing and passing configurations
   to the registration methods `CGAL::pointmatcher::compute_registration_transformation()`
   and `CGAL::pointmatcher::register_point_sets()`.

   \details A configuration corresponds to a component of a configuration module
   of \ref thirdpartylibpointmatcher library. The name and the parameters of any
   configuration for the corresponding registration methods are directly passed
   to \ref thirdpartylibpointmatcher library to be parsed and registered at the
   \ref thirdpartylibpointmatcher side.
   */
struct ICP_config {
  /// The name of the configuration component
  std::string name;

  /// The set of (parameter name, parameter value) pairs as a map
  std::map<std::string, std::string> params;
};

namespace internal {

void dump_invalid_point_matcher_config_exception_msg(const PointMatcherSupport::InvalidElement& err) {
  std::cerr << "ERROR Invalid configuration for PM::ICP, omitting configuration: " << std::endl;
        std::cerr << "   " << err.what() << std::endl;
}

template<typename Scalar, typename NamedParameters1, typename NamedParameters2>
ICP<Scalar>
construct_icp(const NamedParameters1& np1, const NamedParameters2& np2)
{
  typedef PointMatcher<Scalar> PM;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  ICP<Scalar> icp;

  icp.setDefault();

  const ICP_config null_config { "_null_pm_config_in_cgal" };
  const std::vector<ICP_config> null_config_chain { null_config };
  auto is_null_config = [&](const ICP_config& c) { return !c.name.compare(null_config.name); };

  // Config file
  std::istream* pointmatcher_config = choose_parameter(get_parameter(np1, internal_np::pointmatcher_config), nullptr);
  if(pointmatcher_config != nullptr)
  {
    icp.loadFromYaml(*pointmatcher_config);
    return icp;
  }

  // In CGAL, point_set_1 is the reference while point_set_2 is the data
  // However, in pointmatcher, the order is reverse: point_set_1 is the data while point_set_2 is the reference
  // Therefore, filter params from np1 applies to reference data points while params from np2 applies to reading data points

  // np1.point_set_filters -> PM::ReferenceDataPointsFilter
  auto reference_data_points_filter_configs = choose_parameter(get_parameter(np1, internal_np::point_set_filters), null_config_chain);
  if(!reference_data_points_filter_configs.empty() && is_null_config(*reference_data_points_filter_configs.cbegin())) {
    // No config provided: use default values as been set above, do nothing
    ;
  } else {
    // Some config chain is given: clear default values and use the provided configs
    icp.referenceDataPointsFilters.clear();

    for(const auto& conf : reference_data_points_filter_configs)
    {
      try {
        icp.referenceDataPointsFilters.push_back( PM::get().DataPointsFilterRegistrar.create(conf.name, conf.params) );
      } catch(typename PointMatcherSupport::InvalidElement& error) {
        dump_invalid_point_matcher_config_exception_msg(error);
      }
    }
  }
  // np2.point_set_filters -> PM::ReadingDataPointsFilters
  auto reading_data_points_filter_configs = choose_parameter(get_parameter(np2, internal_np::point_set_filters), null_config_chain);
  if(!reading_data_points_filter_configs.empty() && is_null_config(*reading_data_points_filter_configs.cbegin())) {
    // No config provided: use default values as been set above, do nothing
    ;
  } else {
    // Some config chain is given: clear default values and use the provided configs
    icp.readingDataPointsFilters.clear();

    for(const auto& conf : reading_data_points_filter_configs)
    {
      try {
        icp.readingDataPointsFilters.push_back( PM::get().DataPointsFilterRegistrar.create(conf.name, conf.params) );
      } catch(typename PointMatcherSupport::InvalidElement& error) {
        dump_invalid_point_matcher_config_exception_msg(error);
      }
    }
  }

  // Matcher
  auto matcher_config = choose_parameter(get_parameter(np1, internal_np::matcher), null_config);
  if(!is_null_config(matcher_config))
  {
    try {
      icp.matcher = PM::get().MatcherRegistrar.create(matcher_config.name, matcher_config.params);
    } catch(typename PointMatcherSupport::InvalidElement& error) {
      dump_invalid_point_matcher_config_exception_msg(error);
    }
  }

  // Outlier Filters
  auto outlier_filters_config = choose_parameter(get_parameter(np1, internal_np::outlier_filters), null_config_chain);
  if(!outlier_filters_config.empty() && is_null_config(*outlier_filters_config.cbegin())) {
    // No config provided: use default values as been set above, do nothing
    ;
  } else {
    // Some config chain is given: clear default values and use the provided configs
    icp.outlierFilters.clear();

    for(const auto& conf : outlier_filters_config)
    {
      try {
        icp.outlierFilters.push_back( PM::get().OutlierFilterRegistrar.create(conf.name, conf.params) );
      } catch(typename PointMatcherSupport::InvalidElement& error) {
        dump_invalid_point_matcher_config_exception_msg(error);
      }
    }
  }

  // Error Minimizer
  auto error_minimizer_config = choose_parameter(get_parameter(np1, internal_np::error_minimizer), null_config);
  if(!is_null_config(error_minimizer_config))
  {
    try {
      icp.errorMinimizer = PM::get().ErrorMinimizerRegistrar.create(error_minimizer_config.name, error_minimizer_config.params);
    } catch(typename PointMatcherSupport::InvalidElement& error) {
      dump_invalid_point_matcher_config_exception_msg(error);
    }
  }

  // Transformation Checkers
  auto transformation_checkers_config = choose_parameter(get_parameter(np1, internal_np::transformation_checkers), null_config_chain);
  if(!transformation_checkers_config.empty() && is_null_config(*transformation_checkers_config.cbegin())) {
    // No config provided: use default values as been set above, do nothing
    ;
  } else {
    // Some config chain is given: clear default values and use the provided configs
    icp.transformationCheckers.clear();

    for(const auto& conf : transformation_checkers_config)
    {
      try {
        icp.transformationCheckers.push_back( PM::get().TransformationCheckerRegistrar.create(conf.name, conf.params) );
      } catch(typename PointMatcherSupport::InvalidElement& error) {
        dump_invalid_point_matcher_config_exception_msg(error);
      }
    }
  }

  // Inspector
  auto inspector_config = choose_parameter(get_parameter(np1, internal_np::inspector), null_config);
  if(!is_null_config(error_minimizer_config))
  {
    try {
      icp.inspector = PM::get().InspectorRegistrar.create(inspector_config.name, inspector_config.params);
    } catch(typename PointMatcherSupport::InvalidElement& error) {
      dump_invalid_point_matcher_config_exception_msg(error);
    }
  }

  // Logger
  auto logger_config = choose_parameter(get_parameter(np1, internal_np::logger), null_config);
  if(!is_null_config(logger_config))
  {
    try {
      PointMatcherSupport::setLogger( PM::get().LoggerRegistrar.create(logger_config.name, logger_config.params) );
    } catch(typename PointMatcherSupport::InvalidElement& error) {
      dump_invalid_point_matcher_config_exception_msg(error);
    }
  }

  return icp;
}

template<typename Scalar,
         typename PointRange,
         typename PointMap,
         typename VectorMap,
         typename PM_matrix>
void
copy_cgal_points_to_pm_matrix
(const PointRange& prange, PointMap point_map, VectorMap vector_map, PM_matrix& pm_points, PM_matrix& pm_normals)
{
  int idx = 0;
  for(const auto& p : prange)
  {
    // position
    const auto& pos = get(point_map, p);
    pm_points(0, idx) = pos.x();
    pm_points(1, idx) = pos.y();
    pm_points(2, idx) = pos.z();
    pm_points(3, idx) = Scalar(1.);

    // normal
    const auto& normal = get (vector_map, p);
    pm_normals(0, idx) = normal.x();
    pm_normals(1, idx) = normal.y();
    pm_normals(2, idx) = normal.z();

    ++idx;
  }
}

template <class Kernel,
          class PointRange1,
          class PointRange2,
          class PointMap1,
          class PointMap2,
          class VectorMap1,
          class VectorMap2>
std::pair<typename Kernel::Aff_transformation_3, bool>
compute_registration_transformation(const PointRange1& range1, const PointRange2& range2,
                                    PointMap1 point_map1, PointMap2 point_map2,
                                    VectorMap1 vector_map1, VectorMap2 vector_map2,
                                    const typename Kernel::Aff_transformation_3& initial_transform,
                                    ICP<typename Kernel::FT> icp)
{
  using Scalar    = typename Kernel::FT;

  using PM                  = PointMatcher<Scalar>;
  using PM_cloud            = typename PM::DataPoints;
  using PM_matrix           = typename PM::Matrix;
  using PM_transform        = typename PM::Transformation;
  using PM_transform_params = typename PM::TransformationParameters;

  // ref_points: 1, points: 2
  std::size_t nb_ref_points = range1.size();
  std::size_t nb_points     = range2.size();

  PM_matrix ref_points_pos_matrix    = PM_matrix (4, nb_ref_points);
  PM_matrix ref_points_normal_matrix = PM_matrix (3, nb_ref_points);
  PM_matrix points_pos_matrix    = PM_matrix (4, nb_points);
  PM_matrix points_normal_matrix = PM_matrix (3, nb_points);

  // In CGAL, point_set_1 is the reference while point_set_2 is the data

  // convert cgal points to pointmatcher points
  internal::copy_cgal_points_to_pm_matrix<Scalar>(range1,
                                                  point_map1,
                                                  vector_map1,
                                                  ref_points_pos_matrix, // out
                                                  ref_points_normal_matrix); // out

  internal::copy_cgal_points_to_pm_matrix<Scalar>(range2,
                                                  point_map2,
                                                  vector_map2,
                                                  points_pos_matrix, // out
                                                  points_normal_matrix); // out

  auto construct_PM_cloud = [](const PM_matrix& positions, const PM_matrix& normals) -> PM_cloud
  {
    PM_cloud cloud;

    cloud.addFeature("x", positions.row(0));
    cloud.addFeature("y", positions.row(1));
    cloud.addFeature("z", positions.row(2));
    cloud.addFeature("pad", positions.row(3));
    cloud.addDescriptor("normals", normals);

    return cloud;
  };

  PM_cloud ref_cloud = construct_PM_cloud(ref_points_pos_matrix, ref_points_normal_matrix);
  PM_cloud cloud     = construct_PM_cloud(points_pos_matrix,     points_normal_matrix);

  PM_transform_params pm_transform_params = PM_transform_params::Identity(4,4);

  // Convert CGAL transform to pm transform
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      pm_transform_params(i,j) = initial_transform.m(i,j);

  bool converged = false;
  try
  {
                const PM_transform_params prior = pm_transform_params;
                pm_transform_params = icp(cloud, ref_cloud, prior);
    converged = true;
        }
        catch (typename PM::ConvergenceError& error)
        {
                std::cerr << "ERROR CGAL::pointmatcher registration (PM::ICP) failed to converge: " << std::endl;
                std::cerr << "   " << error.what() << std::endl;
    converged = false;
        }

        // Rigid transformation
        std::shared_ptr<PM_transform> transform = PM::get().REG(Transformation).create("RigidTransformation");
  pm_transform_params = transform->correctParameters(pm_transform_params);

  typename Kernel::Aff_transformation_3 cgal_transform
    (pm_transform_params(0,0), pm_transform_params(0,1), pm_transform_params(0,2), pm_transform_params(0,3),
     pm_transform_params(1,0), pm_transform_params(1,1), pm_transform_params(1,2), pm_transform_params(1,3),
     pm_transform_params(2,0), pm_transform_params(2,1), pm_transform_params(2,2), pm_transform_params(2,3));

#ifdef CGAL_POINTMATCHER_VERBOSE
  std::cerr << "Transformation matrix: " << std::endl;
  for (std::size_t i = 0; i < 4; ++ i)
  {
    for (std::size_t j = 0; j < 4; ++ j)
      std::cerr << cgal_transform.coeff(i,j) << " ";
    std::cerr << std::endl;
  }
#endif

  return std::make_pair(cgal_transform, converged);
}

} // end of namespace internal

// point_set_1 is reference while point_set_2 is data
/**
   \ingroup PkgPointSetProcessing3Algorithms

   Computes the registration of `point_set_2` with respect to `point_set_1` and
   returns the corresponding affine transformation.
   Registration is computed using the Iterative Closest Point (ICP) algorithm.

   \note This function requires the \ref thirdpartylibpointmatcher library.

   \tparam PointRange1 is a model of `Range`. The value type of its iterator is
   the key type of the named parameter `point_map` in `NamedParameters1`.
   \tparam PointRange2 is a model of `Range`. The value type of its iterator is
   the key type of the named parameter `point_map` in `NamedParameters2`.

   \param point_set_1 input point range used as reference.
   \param point_set_2 input point range whose registration w.r.t. `point_set_1` will be computed.
   \param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `point_set_1`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange1` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point set `point_set_1`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange1` and whose value type is `geom_traits::Vector_3`}
     \cgalParamNEnd

     \cgalParamNBegin{point_set_filters}
       \cgalParamDescription{a chain of filters to be applied to the point set}
       \cgalParamType{a class model of `Range`. The value type of its iterator must be `ICP_config`.}
       \cgalParamDefault{`RandomSamplingDataPointsFilter`}
       \cgalParamExtra{The chain of filters to be applied to the reference point cloud. The reference
                       point cloud is processed into an intermediate point cloud with the given chain
                       of filters to be used in the alignment procedure. The chain is organized with
                       the forward traversal order of the point set filters range.

                       The chain of point set filters are applied only once at the beginning of the
                       ICP procedure, i.e., before the first iteration of the ICP algorithm.

                       The filters can have several purposes, including but are not limited to
                       i) removal of noisy points which render alignment of point clouds difficult,
                       ii) removal of redundant points so as to speed up alignment, iii) addition
                       of descriptive information to the points such as a surface normal vector,
                       or the direction from the point to the sensor.}
       \cgalParamExtra{Corresponds to `referenceDataPointsFilters` configuration module of \ref thirdpartylibpointmatcher
                       library. The filters should be chosen and set from possible components of
                       the `referenceDataPointsFilters` configuration module.
                       See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
                       for possible configurations.}
     \cgalParamNEnd

     \cgalParamNBegin{matcher}
       \cgalParamDescription{a method used for matching (linking) the points from `point_set_2`,
                             to the points in the reference cloud, `point_set_1`}
       \cgalParamType{a class model of `ICP_config`}
       \cgalParamDefault{`KDTreeMatcher`}
       \cgalParamExtra{Corresponds to the `matcher` configuration module of \ref thirdpartylibpointmatcher
                       library. The matcher should be chosen and set from possible components of
                       the `matcher` configuration module.
                       See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
                       for possible configurations.}
     \cgalParamNEnd

     \cgalParamNBegin{outlier_filters}
       \cgalParamDescription{a chain of filters to be applied to the matched (linked) point clouds after
                             each processing iteration of the ICP algorithm to remove the links which do not
                             correspond to true point correspondences}
       \cgalParamType{a model of `Range`. The value type of its iterator must be `ICP_config`.}
       \cgalParamDefault{`TrimmedDistOutlierFilter`}
       \cgalParamExtra{The outliers are rejected. Points with no link are ignored
                       in the subsequent error minimization step. The chain is organized
                       with the forward traversal order of the outlier filters range.}
       \cgalParamExtra{Corresponds to the `outlierFilters` configuration module of \ref thirdpartylibpointmatcher
                       library. The filters should be chosen and set from possible components of
                       the `outlierFilters` configuration module.
                       See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
                       for possible configurations.}
     \cgalParamNEnd

     \cgalParamNBegin{error_minimizer}
       \cgalParamDescription{an error minimizer that computes a transformation matrix such as to minimize
                             the error between the point sets}
       \cgalParamType{a class model of `ICP_config`}
       \cgalParamDefault{`PointToPlaneErrorMinimizer`}
       \cgalParamExtra{Corresponds to the `errorMinimizer` configuration module of \ref thirdpartylibpointmatcher
                       library. The error minimizer should be chosen and set from possible components of
                       the `errorMinimizer` configuration module.
                       See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
                       for possible configurations.}
     \cgalParamNEnd

     \cgalParamNBegin{transformation_checkers}
       \cgalParamDescription{a chain of transformation checkers}
       \cgalParamType{a class model of `Range`. The value type of its iterator must be `ICP_config`.}
       \cgalParamDefault{`CounterTransformationChecker` and `DifferentialTransformationChecker`}
       \cgalParamExtra{The chain is organized with the forward traversal order of the transformation checkers range.}
       \cgalParamExtra{A transformation checker can stop the iteration depending on the conditions it defines.}
       \cgalParamExtra{Corresponds to the `transformationCheckers` configuration module of \ref thirdpartylibpointmatcher
                       library. The transformation checkers should be chosen and set from possible components of
                       the `transformationCheckers` configuration module.
                       See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
                       for possible configurations.}
     \cgalParamNEnd

     \cgalParamNBegin{inspector}
       \cgalParamDescription{an inspector that enables logging data at different steps for analysis.}
       \cgalParamType{a class model of `ICP_config`}
       \cgalParamDefault{`NullInspector`}
       \cgalParamExtra{Inspectors typically provide deeper scrutiny than the logger.}
       \cgalParamExtra{Corresponds to the `inspector` configuration module of \ref thirdpartylibpointmatcher
                       library. The inspector should be chosen and set from possible components of
                       the `inspector` configuration module.
                       See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
                       for possible configurations.}
     \cgalParamNEnd

     \cgalParamNBegin{logger}
       \cgalParamDescription{a method for logging information regarding the registration process
                             outputted by \ref thirdpartylibpointmatcher library}
       \cgalParamType{a class model of `ICP_config`}
       \cgalParamDefault{`NullLogger`}
       \cgalParamExtra{The logs generated by CGAL library does not get effected by this configuration.}
       \cgalParamExtra{Corresponds to the `logger` configuration module of \ref thirdpartylibpointmatcher
                       library. The logger should be chosen and set from possible components of
                       the `logger` configuration module.
                       See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
                       for possible configurations.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `point_set_2`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange2` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point set `point_set_2`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange2` and whose value type is `geom_traits::Vector_3`}
     \cgalParamNEnd

     \cgalParamNBegin{point_set_filters}
       \cgalParamDescription{a chain of filters to be applied to the point set}
       \cgalParamType{a class model of `Range`. The value type of its iterator must be `ICP_config`.}
       \cgalParamDefault{`SamplingSurfaceNormalDataPointsFilter`}
       \cgalParamExtra{The chain of filters to be applied to the point cloud `point_set_2`. The
                       point cloud is processed into an intermediate point cloud with the given chain
                       of filters to be used in the alignment procedure. The chain is organized with
                       the forward traversal order of the point set filters range.

                       The chain of point set filters are applied only once at the beginning of the
                       ICP procedure, i.e., before the first iteration of the ICP algorithm.

                       The filters can have several purposes, including but are not limited to
                       i) removal of noisy points which render alignment of point clouds difficult,
                       ii) removal of redundant points so as to speed up alignment, iii) addition
                       of descriptive information to the points such as a surface normal vector,
                       or the direction from the point to the sensor.}
       \cgalParamExtra{Corresponds to the `readingDataPointsFilters` configuration module of \ref thirdpartylibpointmatcher
                       library. The filters should be chosen and set from possible components of
                       the `readingDataPointsFilters` configuration module.
                       See <a href="https://libpointmatcher.readthedocs.io/en/latest/Configuration/#configuration-of-an-icp-chain">libpointmatcher documentation</a>
                       for possible configurations.}
     \cgalParamNEnd

     \cgalParamNBegin{transformation}
       \cgalParamDescription{an affine transformation that is used as the initial transformation for `point_set_2`}
       \cgalParamType{`CGAL::Aff_transformation_3`}
       \cgalParamDefault{the identity transformation}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \return a pair containing the affine transformation that should be applied
   to `point_set_2` to make it registered w.r.t. `point_set_1` and the
   boolean value indicating if the registration converged. The second
   of the pair is `true` if converged, `false` otherwise. A log why it failed to
   converge is written to `std::cerr` if the registration cannot converge.
*/
template <class PointRange1, class PointRange2,
          class NamedParameters1, class NamedParameters2>
#ifdef DOXYGEN_RUNNING
std::pair<geom_traits::Aff_transformation_3, bool>
#else
std::pair<typename CGAL::Point_set_processing_3::GetK<PointRange1, NamedParameters1>
  ::Kernel::Aff_transformation_3, bool>
#endif
compute_registration_transformation (const PointRange1& point_set_1, const PointRange2& point_set_2,
                                     const NamedParameters1& np1, const NamedParameters2& np2)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  namespace PSP = CGAL::Point_set_processing_3;

  // property map types
  typedef typename CGAL::GetPointMap<PointRange1, NamedParameters1>::type PointMap1;
  typedef typename CGAL::GetPointMap<PointRange2, NamedParameters2>::type PointMap2;
  CGAL_static_assertion_msg((boost::is_same< typename boost::property_traits<PointMap1>::value_type,
                                             typename boost::property_traits<PointMap2>::value_type> ::value),
                            "The point type of input ranges must be the same");

  typedef typename PSP::GetNormalMap<PointRange1, NamedParameters1>::type NormalMap1;
  typedef typename PSP::GetNormalMap<PointRange2, NamedParameters2>::type NormalMap2;
  CGAL_static_assertion_msg((boost::is_same< typename boost::property_traits<NormalMap1>::value_type,
                                             typename boost::property_traits<NormalMap2>::value_type> ::value),
                            "The vector type of input ranges must be the same");

  typedef typename PSP::GetK<PointRange1, NamedParameters1>::Kernel Kernel;
  typedef typename Kernel::FT Scalar;
  typedef typename Kernel::Aff_transformation_3 Transformation;

  PointMap1 point_map1 = choose_parameter(get_parameter(np1, internal_np::point_map), PointMap1());
  NormalMap1 normal_map1 = choose_parameter(get_parameter(np1, internal_np::normal_map), NormalMap1());
  PointMap2 point_map2 = choose_parameter(get_parameter(np2, internal_np::point_map), PointMap2());
  NormalMap2 normal_map2 = choose_parameter(get_parameter(np2, internal_np::normal_map), NormalMap2());

  // initial transformation
  Transformation initial_transformation
    = choose_parameter(get_parameter(np2, internal_np::transformation), Transformation(Identity_transformation()));

  return internal::compute_registration_transformation<Kernel>(point_set_1, point_set_2,
                                                               point_map1, point_map2,
                                                               normal_map1, normal_map2,
                                                               initial_transformation,
                                                               internal::construct_icp<Scalar>(np1, np2));
}

// convenience overloads
template <class PointRange1, class PointRange2,
          class NamedParameters1>
std::pair<typename CGAL::Point_set_processing_3::GetK<PointRange1, NamedParameters1>
  ::Kernel::Aff_transformation_3, bool>
compute_registration_transformation(const PointRange1& point_set_1, const PointRange2& point_set_2,
      const NamedParameters1& np1)
{
  namespace params = CGAL::Point_set_processing_3::parameters;
  return compute_registration_transformation(point_set_1, point_set_2, np1, params::all_default(point_set_1));
}

template <class PointRange1, class PointRange2>
std::pair<typename CGAL::Point_set_processing_3::GetK<PointRange1,
          Named_function_parameters<bool, internal_np::all_default_t> >
  ::Kernel::Aff_transformation_3, bool>
compute_registration_transformation(const PointRange1& point_set_1, const PointRange2& point_set_2)
{
  namespace params = CGAL::Point_set_processing_3::parameters;
  return compute_registration_transformation(point_set_1, point_set_2,
                                             params::all_default(point_set_1),
                                             params::all_default(point_set_2));
}


} } // end of namespace CGAL::pointmatcher

#endif // CGAL_LINKED_WITH_POINTMATCHER

#endif // CGAL_POINTMATCHER_COMPUTE_REGISTRATION_TRANSFORMATION_H
