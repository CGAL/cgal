// Copyright (c) 2019  GeometryFactory(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Sebastien Loriot, Necip Fazil Yildiran

#ifndef CGAL_OPENGR_COMPUTE_REGISTRATION_TRANSFORMATION_H
#define CGAL_OPENGR_COMPUTE_REGISTRATION_TRANSFORMATION_H

#include <CGAL/license/Point_set_processing_3.h>

#if defined(CGAL_LINKED_WITH_OPENGR) || defined(DOXYGEN_RUNNING)

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Iterator_range.h>

#include <boost/type_traits/is_same.hpp>
#include <boost/range/iterator_range.hpp>

#include <gr/algorithms/FunctorSuper4pcs.h>
#include <gr/algorithms/PointPairFilter.h>

#include <Eigen/Dense>

namespace CGAL {

namespace OpenGR {

template<typename Kernel>
using Options = typename gr::Match4pcsBase<gr::FunctorSuper4PCS,
                                           gr::Point3D<typename Kernel::FT>,
                                           gr::DummyTransformVisitor,
                                           gr::AdaptivePointFilter,
                                           gr::AdaptivePointFilter::Options>::OptionsType;

namespace internal {

template <typename Scalar, typename InputRange, typename PointMap, typename VectorMap>
struct CGAL_range_and_pmaps_to_opengr_point3d_range
{
  typedef typename InputRange::const_iterator::value_type argument_type;
  typedef gr::Point3D<Scalar> result_type;
  typedef typename result_type::VectorType vector_type;

  PointMap point_map;
  VectorMap normal_map;

  CGAL_range_and_pmaps_to_opengr_point3d_range (PointMap point_map, VectorMap normal_map)
    : point_map (point_map), normal_map (normal_map)
  { }

  result_type operator() (const argument_type& arg) const
  {
    const auto& p = get (point_map, arg);
    const auto& n = get (normal_map, arg);

    result_type out (p.x(), p.y(), p.z());
    out.set_normal ( vector_type(n.x(), n.y(), n.z()) );

    return out;
  }
};

template <class Kernel,
          class PointRange1,
          class PointRange2,
          class PointMap1,
          class PointMap2,
          class VectorMap1,
          class VectorMap2>
std::pair<typename Kernel::Aff_transformation_3, double>
compute_registration_transformation(const PointRange1& range1,    const PointRange2& range2,
                                    PointMap1 point_map1,   PointMap2 point_map2,
                                    VectorMap1 vector_map1, VectorMap2 vector_map2,
                                    Options<Kernel>& options)
{
  typedef gr::Point3D<typename Kernel::FT> PointType;

  namespace GR=gr;

  // TODO: see if should allow user to change those types
  typedef Eigen::Matrix<typename PointType::Scalar, 4, 4> MatrixType;
  typedef gr::UniformDistSampler<PointType> SamplerType;
  typedef gr::DummyTransformVisitor TrVisitorType;
  typedef gr::Match4pcsBase<gr::FunctorSuper4PCS,
                            PointType,
                            TrVisitorType,
                            gr::AdaptivePointFilter,
                            gr::AdaptivePointFilter::Options> MatcherType;

  MatrixType mat (MatrixType::Identity());
  SamplerType sampler;
  TrVisitorType visitor;

  // Unary functions that convert value_type of point ranges to gr::Point3D
  CGAL_range_and_pmaps_to_opengr_point3d_range<typename Kernel::FT, PointRange1, PointMap1, VectorMap1> // TODO: remove deductible ones
    unary_function_1 (point_map1, vector_map1);

  CGAL_range_and_pmaps_to_opengr_point3d_range<typename Kernel::FT, PointRange2, PointMap2, VectorMap2> // TODO: remove deductible ones
    unary_function_2 (point_map2, vector_map2);

  auto gr_point_range_1 = boost::make_iterator_range(
    boost::make_transform_iterator (range1.begin(), unary_function_1),
    boost::make_transform_iterator (range1.end(),   unary_function_1));

  auto gr_point_range_2 = boost::make_iterator_range(
    boost::make_transform_iterator (range2.begin(), unary_function_2),
    boost::make_transform_iterator (range2.end(),   unary_function_2));

  // logger
  GR::Utils::Logger logger(GR::Utils::NoLog);

  // matcher
  MatcherType matcher(options, logger);
  double score =
    matcher.ComputeTransformation(gr_point_range_1, gr_point_range_2, mat, sampler, visitor);

#ifdef CGAL_OPENGR_VERBOSE
  std::cerr << "Transformation matrix: " << std::endl;
  for (std::size_t i = 0; i < 4; ++ i)
  {
    for (std::size_t j = 0; j < 4; ++ j)
      std::cerr << mat.coeff(i,j) << " ";
    std::cerr << std::endl;
  }
#endif

  typename Kernel::Aff_transformation_3 cgal_trsf(
    mat.coeff(0,0), mat.coeff(0,1), mat.coeff(0,2), mat.coeff(0,3),
    mat.coeff(1,0), mat.coeff(1,1), mat.coeff(1,2), mat.coeff(1,3),
    mat.coeff(2,0), mat.coeff(2,1), mat.coeff(2,2), mat.coeff(2,3));

  return std::make_pair(cgal_trsf, score);
}

} // end of namespace internal

/**
   \ingroup PkgPointSetProcessing3Algorithms

   Computes the registration of `point_set_2` with respect to `point_set_1` and
   returns the corresponding affine transformation along with the registration
   score.

   Registration is computed using the Super4PCS algorithm \cgalCite{cgal:mam-sffgp-14}.

   \note This function requires the \ref thirdpartyOpenGR library.

   \warning Although this may seem counter-intuitive, if one of the
   two point set matches only a small section of the other one, it is
   advised to _use the small point set as reference_ instead of the
   big one. The reason is that the reference point set is used to
   construct a base that is sought after in the other point set: if
   the big point set is used as reference, chances are the constructed
   base will not be present in the small point set.

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
       \cgalParamDefault{Normals are computed and stored internally.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd

     \cgalParamNBegin{number_of_samples}
       \cgalParamDescription{size of the subset of input points used to compute registration}
       \cgalParamType{unsigned int}
       \cgalParamDefault{`200`}
       \cgalParamExtra{Input clouds are sub-sampled prior exploration, to ensure fast computations.
                       Super4PCS has a linear complexity w.r.t. the number of input samples,
                       allowing to use larger values than 4PCS. Simple geometry with large overlap
                       can be matched with only 200 samples. However, with Super4PCS, smaller details
                       can be used during the process by using up to thousands of points.
                       There is no theoretical limit to this parameter, however using too large values
                       leads to very a large congruent set, which requires more time and memory
                       to be explored. Using a large number of samples is recommended when:
                       geometrical details are required to perform the matching, for instance
                       to disambiguate between several similar configurations; the clouds
                       have a very low overlap: using a too sparse sampling can prevent
                       to have samples in the overlapping area, causing the algorithm to fail;
                       the clouds are very noisy, and require a dense sampling.
                       Note that Super4PCS is a global registration algorithm, which
                       finds a good approximate of the rigid transformation aligning too
                       clouds. Increasing the number of samples in order to get a fine
                       registration is not optimal: it is usually faster to use less samples, and
                       refine the transformation using a local algorithm, like the ICP, or its
                       variant SparseICP.}
     \cgalParamNEnd

     \cgalParamNBegin{maximum_normal_deviation}
       \cgalParamDescription{angle threshold (in degrees) used to filter pairs of points
                             according to their normal consistency}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{`90.00`}
       \cgalParamExtra{Small values decrease computation time but may also decrease the quality
                       if pairs of points that should match have a normal deviation higher than the threshold.}
     \cgalParamNEnd

     \cgalParamNBegin{accuracy}
       \cgalParamDescription{registration accuracy (delta in the paper)}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{`5.00`}
       \cgalParamExtra{Setting a small value means that the two clouds needs to be very
                       close to be considered as well aligned. It is expressed in scene units. A
                       simple way to understand its impact is to consider the computation of the
                       Largest Common Pointset (LCP), the metric used to verify how much the
                       clouds are aligned. For each transformation matrix produced by Super4PCS,
                       we compute the LCP measure by considering a shell around the reference
                       cloud, and count the percentage of points of the target cloud lying in the
                       shell. The thickness of the shell is defined by the parameter
                       delta.}
     \cgalParamNEnd

     \cgalParamNBegin{overlap}
       \cgalParamDescription{ratio of expected overlap between the two point sets:
                             it is ranging between `0` (no overlap) to `1` (100% overlap)}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{`0.20`}
       \cgalParamExtra{The overlap parameter controls the size of the basis used for
                       registration. Usually, the larger the overlap, the faster the
                       algorithm. When the overlap is unknown, a simple way to set this parameter
                       is to start from 100% overlap, and decrease the value until obtaining a
                       good result. Using too small values will slow down the algorithm, and
                       reduce the accuracy of the result.}
     \cgalParamNEnd

     \cgalParamNBegin{maximum_running_time}
       \cgalParamDescription{maximum number of seconds after which the algorithm terminates.}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{`1000`}
       \cgalParamExtra{Super4PCS explores the transformation space to align the two input clouds.
                       Since the exploration is performed randomly, it is recommended
                       to use a large time value to explore the whole space.}
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
       \cgalParamDescription{a property map associating normals to the elements of the poing set `point_set_2`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange2` and whose value type is `geom_traits::Vector_3`}
       \cgalParamDefault{Normals are computed and stored internally.}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \return a pair containing the affine transformation that should be applied
   to `point_set_2` to make it registered w.r.t. `point_set_1` and the
   registration score.
*/
template <class PointRange1, class PointRange2,
          class NamedParameters1, class NamedParameters2>
#ifdef DOXYGEN_RUNNING
std::pair<geom_traits::Aff_transformation_3, double>
#else
std::pair<typename CGAL::Point_set_processing_3::GetK<PointRange1, NamedParameters1>
  ::Kernel::Aff_transformation_3, double>
#endif
compute_registration_transformation (const PointRange1& point_set_1, const PointRange2& point_set_2,
                                     const NamedParameters1& np1, const NamedParameters2& np2)
{
  namespace PSP = CGAL::Point_set_processing_3;
  namespace GR = gr;
  using parameters::choose_parameter;
  using parameters::get_parameter;

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

  PointMap1 point_map1 = choose_parameter(get_parameter(np1, internal_np::point_map), PointMap1());
  NormalMap1 normal_map1 = choose_parameter(get_parameter(np1, internal_np::normal_map), NormalMap1());
  PointMap2 point_map2 = choose_parameter(get_parameter(np2, internal_np::point_map), PointMap2());
  NormalMap2 normal_map2 = choose_parameter(get_parameter(np2, internal_np::normal_map), NormalMap2());

  Options<Kernel> options;
  options.sample_size = choose_parameter(get_parameter(np1, internal_np::number_of_samples), 200);
  options.delta = choose_parameter(get_parameter(np1, internal_np::accuracy), 5.00);
  options.max_time_seconds = choose_parameter(get_parameter(np1, internal_np::maximum_running_time), 1000);
  options.max_normal_difference = choose_parameter(get_parameter(np1, internal_np::maximum_normal_deviation), 90.);

  bool overlap_ok = options.configureOverlap (choose_parameter(get_parameter(np1, internal_np::overlap), 0.20));
  CGAL_USE (overlap_ok);
  CGAL_assertion_msg (overlap_ok, "Invalid overlap configuration.");

  return internal::compute_registration_transformation<Kernel>(point_set_1, point_set_2,
                                                               point_map1, point_map2,
                                                               normal_map1, normal_map2,
                                                               options);

}

// convenience overloads
template <class PointRange1, class PointRange2,
          class NamedParameters1>
std::pair<typename CGAL::Point_set_processing_3::GetK<PointRange1, NamedParameters1>
  ::Kernel::Aff_transformation_3, double>
compute_registration_transformation(const PointRange1& point_set_1, PointRange2& point_set_2,
      const NamedParameters1& np1)
{
  namespace params = CGAL::Point_set_processing_3::parameters;
  return compute_registration_transformation(point_set_1, point_set_2, np1, params::all_default(point_set_1));
}

template <class PointRange1, class PointRange2>
std::pair<typename CGAL::Point_set_processing_3::GetK<PointRange1,
          Named_function_parameters<bool, internal_np::all_default_t> >
  ::Kernel::Aff_transformation_3, double>
compute_registration_transformation(const PointRange1& point_set_1, PointRange2& point_set_2)
{
  namespace params = CGAL::Point_set_processing_3::parameters;
  return compute_registration_transformation(point_set_1, point_set_2,
                                             params::all_default(point_set_1),
                                             params::all_default(point_set_2));
}

} } // end of namespace CGAL::OpenGR

#endif // CGAL_LINKED_WITH_OPENGR

#endif // CGAL_OPENGR_COMPUTE_REGISTRATION_TRANSFORMATION_H
