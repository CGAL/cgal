// Copyright (c) 2018  GeometryFactory(France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Sebastien Loriot

#ifndef CGAL_OPENGR_COMPUTE_REGISTRATION_TRANSFORMATION_H
#define CGAL_OPENGR_COMPUTE_REGISTRATION_TRANSFORMATION_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/type_traits/is_same.hpp>

#include <gr/algorithms/FunctorSuper4pcs.h>
#include <gr/algorithms/PointPairFilter.h>

#include <Eigen/Dense>

namespace CGAL {

namespace OpenGR {

typedef gr::Match4pcsBase<gr::FunctorSuper4PCS,
                          gr::DummyTransformVisitor,
                          gr::AdaptivePointFilter,
                          gr::AdaptivePointFilter::Options>::OptionsType Options;

namespace internal {

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
                                    Options& options)
{
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  namespace GR=gr;

  std::vector<GR::Point3D> set1, set2;
  std::vector<GR::Point3D::VectorType> normals1, normals2;

  // TODO: see if should allow user to change those types
  typedef Eigen::Matrix<gr::Point3D::Scalar, 4, 4> MatrixType;
  typedef gr::UniformDistSampler SamplerType;
  typedef gr::DummyTransformVisitor TrVisitorType;
  typedef gr::Match4pcsBase<gr::FunctorSuper4PCS,
                            TrVisitorType,
                            gr::AdaptivePointFilter,
                            gr::AdaptivePointFilter::Options> MatcherType;

  MatrixType mat (MatrixType::Identity());
  SamplerType sampler;
  TrVisitorType visitor;

  // copy points and normal
  const std::size_t nbpt1 = range1.size();
  set1.reserve(nbpt1);
  normals1.reserve(nbpt1);
  for (typename PointRange1::const_iterator it=range1.begin(),
                                            end=range1.end(); it!=end; ++it)
  {
    const Point_3& p = get(point_map1, *it);
    const Vector_3& v = get(vector_map1, *it);
    set1.push_back(GR::Point3D(p.x(), p.y(), p.z()));
    normals1.push_back(GR::Point3D::VectorType(v.x(), v.y(), v.z()));
  }

  const std::size_t nbpt2 = range2.size();
  set2.reserve(nbpt2);
  normals2.reserve(nbpt2);
  for (typename PointRange1::const_iterator it=range2.begin(),
                                            end=range2.end(); it!=end; ++it)
  {
    const Point_3& p = get(point_map2, *it);
    const Vector_3& v = get(vector_map2, *it);
    set2.push_back(GR::Point3D(p.x(), p.y(), p.z()));
    normals2.push_back(GR::Point3D::VectorType(v.x(), v.y(), v.z()));
  }

  // logger
  GR::Utils::Logger logger(GR::Utils::NoLog);

  // matcher
  MatcherType matcher(options, logger);
  double score =
    matcher.ComputeTransformation(set1, set2, mat, sampler, visitor);

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

     \cgalParamBegin{number_of_samples} size of the subset of input points used
     to compute registration. Input clouds are sub-sampled prior exploration,
     to ensure fast computations. Super4PCS has a linear complexity w.r.t. the
     number of input samples, allowing to use larger values than 4PCS. Simple
     geometry with large overlap can be matched with only 200 samples. However,
     with Super4PCS, smaller details can be used during the process by using up
     to thousands of points. There is no theoretical limit to this parameter,
     however using too large values leads to very a large congruent set, which
     requires more time and memory to be explored. Using a large number of
     samples is recommended when: geometrical details are required to perform
     the matching, for instance to disambiguate between several similar
     configurations; the clouds have a very low overlap: using a too sparse
     sampling can prevent to have samples in the overlapping area, causing the
     algorithm to fail; the clouds are very noisy, and require a dense
     sampling. Note that Super4PCS is a global registration algorithm, which
     finds a good approximate of the rigid transformation aligning too
     clouds. Increasing the number of samples in order to get a fine
     registration is not optimal: it is usually faster to use less samples, and
     refine the transformation using a local algorithm, like the ICP, or its
     variant SparseICP.\cgalParamEnd

     \cgalParamBegin{accuracy} registration accuracy (delta in the
     paper). Setting a small value means that the two clouds needs to be very
     close to be considered as well aligned. It is expressed in scene units. A
     simple way to understand its impact is to consider the computation of the
     Largest Common Pointset (LCP), the metric used to verify how much the
     clouds are aligned. For each transformation matrix produced by Super4PCS,
     we compute the LCP measure by considering a shell around the reference
     cloud, and count the % of points of the target cloud lying in the
     shell. The thickness of the shell is defined by the parameter
     delta.\cgalParamEnd

     \cgalParamBegin{overlap} ratio of expected overlap between the two point
     sets: it is ranging between 0 (no overlap) to 1 (100% overlap). The
     overlap parameter controls the size of the basis used for
     registration. Usually, the larger the overlap, the faster the
     algorithm. When the overlap is unknown, a simple way to set this parameter
     is to start from 100% overlap, and decrease the value until obtaining a
     good result. Using too small values will slow down the algorithm, and
     reduce the accuracy of the result.\cgalParamEnd

     \cgalParamBegin{maximum_running_time} maximum number of seconds after
     which the algorithm stops. Super4PCS explores the transformation space to
     align the two input clouds. Since the exploration is performed randomly,
     it is recommended to use a large time value to explore the whole space
     (e.g., 1000).\cgalParamEnd

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
  using boost::choose_param;
  using boost::get_param;

  // property map types
  typedef typename PSP::GetPointMap<PointRange1, NamedParameters1>::type PointMap1;
  typedef typename PSP::GetPointMap<PointRange2, NamedParameters2>::type PointMap2;
  CGAL_static_assertion_msg((boost::is_same< typename boost::property_traits<PointMap1>::value_type,
                                             typename boost::property_traits<PointMap2>::value_type> ::value),
                            "The point type of input ranges must be the same");

  typedef typename PSP::GetNormalMap<PointRange1, NamedParameters1>::type NormalMap1;
  typedef typename PSP::GetNormalMap<PointRange2, NamedParameters2>::type NormalMap2;
  CGAL_static_assertion_msg((boost::is_same< typename boost::property_traits<NormalMap1>::value_type,
                                             typename boost::property_traits<NormalMap2>::value_type> ::value),
                            "The vector type of input ranges must be the same");

  typedef typename PSP::GetK<PointRange1, NamedParameters1>::Kernel Kernel;

  PointMap1 point_map1 = choose_param(get_param(np1, internal_np::point_map), PointMap1());
  NormalMap1 normal_map1 = choose_param(get_param(np1, internal_np::normal_map), NormalMap1());
  PointMap1 point_map2 = choose_param(get_param(np2, internal_np::point_map), PointMap2());
  NormalMap2 normal_map2 = choose_param(get_param(np2, internal_np::normal_map), NormalMap2());

  Options options;
  options.sample_size = choose_param(get_param(np1, internal_np::number_of_samples), 200);
  options.delta = choose_param(get_param(np1, internal_np::accuracy), 5.00);
  options.max_time_seconds = choose_param(get_param(np1, internal_np::maximum_running_time), 1000);
  bool overlap_ok = options.configureOverlap (choose_param(get_param(np1, internal_np::overlap), 0.20));
  CGAL_USE (overlap_ok);
  CGAL_static_assertion_msg (overlap_ok, "Invalid overlap configuration.");

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
          cgal_bgl_named_params<bool, internal_np::all_default_t> >
  ::Kernel::Aff_transformation_3, double>
compute_registration_transformation(const PointRange1& point_set_1, PointRange2& point_set_2)
{
  namespace params = CGAL::Point_set_processing_3::parameters;
  return compute_registration_transformation(point_set_1, point_set_2,
                                             params::all_default(point_set_1),
                                             params::all_default(point_set_2));
}

} } // end of namespace CGAL::OpenGR

#endif // CGAL_OPENGR_COMPUTE_REGISTRATION_TRANSFORMATION_H
