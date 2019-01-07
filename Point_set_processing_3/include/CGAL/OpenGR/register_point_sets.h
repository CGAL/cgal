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

#ifndef CGAL_OPENGR_REGISTER_POINT_SETS_H
#define CGAL_OPENGR_REGISTER_POINT_SETS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/OpenGR/compute_registration_transformation.h>

#include <boost/type_traits/is_same.hpp>

#include <Eigen/Dense>

namespace CGAL {

namespace OpenGR {

namespace internal {

template <class Kernel,
          class PointRange1,
          class PointRange2,
          class PointMap1,
          class PointMap2,
          class VectorMap1,
          class VectorMap2>
double
register_point_sets(const PointRange1& range1,    PointRange2& range2,
                    PointMap1 point_map1,   PointMap2 point_map2,
                    VectorMap1 vector_map1, VectorMap2 vector_map2,
                    Options& options)
{
  std::pair<typename Kernel::Aff_transformation_3, double> res =
    compute_registration_transformation<Kernel>(range1, range2,
                                                point_map1, point_map2,
                                                vector_map1, vector_map2,
                                                options);

  // update CGAL points
  for (typename PointRange2::iterator it=range2.begin(),
                                      end=range2.end(); it!=end; ++it)
  {
    put(point_map2, *it, get(point_map2, *it).transform(res.first));
  }

  return res.second;
}

} // end of namespace internal

/**
   \ingroup PkgPointSetProcessing3Algorithms
   
   Computes the registration of `point_set_2` with respect to `point_set_1` and
   applies it.

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
     to compute registration.\cgalParamEnd
     \cgalParamBegin{accuracy} registration accuracy expressed in scene
     units.\cgalParamEnd
     \cgalParamBegin{overlap} ratio of expected overlap between the two point
     sets (between 0 for no overlap to 1 for a full overlap).\cgalParamEnd
     \cgalParamBegin{maximum_running_time} number of seconds after which the
     algorithm is forced to stop and returns the best solution found so
     far.\cgalParamEnd
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

   \return the registration score.
*/
template <class PointRange1, class PointRange2,
          class NamedParameters1, class NamedParameters2>
double
register_point_sets (const PointRange1& point_set_1, PointRange2& point_set_2,
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

  return internal::register_point_sets<Kernel>(point_set_1, point_set_2,
                                               point_map1, point_map2,
                                               normal_map1, normal_map2,
                                               options);
}

// convenience overloads
template <class PointRange1, class PointRange2,
          class NamedParameters1>
double
register_point_sets(const PointRange1& point_set_1, PointRange2& point_set_2,
                    const NamedParameters1& np1)
{
  namespace params = CGAL::Point_set_processing_3::parameters;
  return register_point_sets(point_set_1, point_set_2, np1, params::all_default(point_set_1));
}

template <class PointRange1, class PointRange2>
double
register_point_sets(const PointRange1& point_set_1, PointRange2& point_set_2)
{
  namespace params = CGAL::Point_set_processing_3::parameters;
  return register_point_sets(point_set_1, point_set_2,
                             params::all_default(point_set_1),
                             params::all_default(point_set_2));
}

} } // end of namespace CGAL::OpenGR

#endif // CGAL_OPENGR_REGISTER_POINT_SETS_H
