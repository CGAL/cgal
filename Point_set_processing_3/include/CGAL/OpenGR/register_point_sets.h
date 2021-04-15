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

#ifndef CGAL_OPENGR_REGISTER_POINT_SETS_H
#define CGAL_OPENGR_REGISTER_POINT_SETS_H

#include <CGAL/license/Point_set_processing_3.h>

#if defined(CGAL_LINKED_WITH_OPENGR) || defined(DOXYGEN_RUNNING)

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
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
                    Options<Kernel>& options)
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

   Registration is computed using the Super4PCS algorithm
   \cgalCite{cgal:mam-sffgp-14}. Parameters documentation is copy-pasted from [the official documentation of OpenGR](https://storm-irit.github.io/OpenGR/a00012.html). For more details on this method, please refer to it.

   \note This function requires the \ref thirdpartyOpenGR library.

   \warning Although this may seem counter-intuitive, if one of the two
   point set matches only a small section of the other one, it is
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
  PointMap1 point_map2 = choose_parameter(get_parameter(np2, internal_np::point_map), PointMap2());
  NormalMap2 normal_map2 = choose_parameter(get_parameter(np2, internal_np::normal_map), NormalMap2());

  Options<Kernel> options;
  options.sample_size = choose_parameter(get_parameter(np1, internal_np::number_of_samples), 200);
  options.delta = choose_parameter(get_parameter(np1, internal_np::accuracy), 5.00);
  options.max_time_seconds = choose_parameter(get_parameter(np1, internal_np::maximum_running_time), 1000);
  options.max_normal_difference = choose_parameter(get_parameter(np1, internal_np::maximum_normal_deviation), 90.);

  bool overlap_ok = options.configureOverlap (choose_parameter(get_parameter(np1, internal_np::overlap), 0.20));
  CGAL_USE (overlap_ok);
  CGAL_assertion_msg (overlap_ok, "Invalid overlap configuration.");

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

#endif // CGAL_LINKED_WITH_OPENGR

#endif // CGAL_OPENGR_REGISTER_POINT_SETS_H
