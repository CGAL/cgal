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

/*!
 * TODO: document me
 */
template <class PointRange1, class PointRange2,
          class NamedParameters1, class NamedParameters2>
std::pair<typename CGAL::Point_set_processing_3::GetK<PointRange1, NamedParameters1>
  ::Kernel::Aff_transformation_3, double>
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

  Options options = choose_param(get_param(np1, internal_np::opengr_options), Options());

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
