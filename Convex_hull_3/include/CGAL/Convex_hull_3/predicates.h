// Copyright (c) 2022  INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Samuel Hornus and Sébastien Loriot
//

#ifndef CGAL_CONVEX_HULL_3_PREDICATES_H
#define CGAL_CONVEX_HULL_3_PREDICATES_H

namespace CGAL {

namespace Convex_hull_3 {

#include <CGAL/Named_function_parameters.h>

/**
* \ingroup PkgConvexHull3Predicates
*
* indicates if the convex hulls of two point sets intersect or not.
*
* @tparam PointRange1 a model of the concept `RandomAccessContainer`
* whose value type is a point type from a CGAL %Kernel
* @tparam PointRange2 a model of the concept `RandomAccessContainer`
* whose value type is the point type from a CGAL %Kernel
* @tparam NamedParameters_1 a sequence of \ref bgl_namedparameters "Named Parameters"
* @tparam NamedParameters_2 a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param r1 first range of points whose convex hull is considered in the do-intersect test
* @param r2 first range of points whose convex hull is considered in the do-intersect test
* @param np_1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
* @param np_2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of the range `r1` (`r2)}
*     \cgalParamType{a model of `ReadablePropertyMap` whose value types are the same for `np1` and `np2`}
*     \cgalParamDefault{`CGAL::Identity_property_map`}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{`np1` only}
*   \cgalParamNEnd
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{if not `0` (no limit), indicates the maximum number of iterations that the algorithm is allowed to do.
*                           If this value is not `0`, then an intersection might be reported even if the convex hulls does not intersect.
                            However, if the convex hulls are reported not to intersect, this is guaranteed.}
*     \cgalParamType{an positive integer convertible to `std::size_t`}
*     \cgalParamDefault{`0`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
*/
template <class PointRange1, class PointRange2,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
bool do_intersect(const PointRange1& r1, const PointRange2& r2,
                  const NamedParameters1 np1 = parameters::default_values(),
                  const NamedParameters2 np2 = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef Point_set_processing_3_np_helper<ReferencePointRange, NamedParameters1> NP_helper1;
  typedef typename NP_helper1::Const_point_map PointMap1;
  typedef typename NP_helper1::Geom_traits Geom_traits;f
  
  typedef Point_set_processing_3_np_helper<ReferencePointRange, NamedParameters2> NP_helper2;
  typedef typename NP_helper2::Const_point_map PointMap2;

  typedef typename boost::property_traits<PointMap1>::value_type Point_3;
  CGAL_static_assertion((std::is_same<Point_3, typename boost::property_traits<PointMap2>::value_type>::value));

  PointMap1 point_map1 = NP_helper1::get_const_point_map(ref_points, np1);
  PointMap2 point_map2 = NP_helper2::get_const_point_map(points, np2);
}

template <class PointRange1, class PointRange2,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
FT
#else
typename Point_set_processing_3_np_helper<ReferencePointRange, NamedParameters1>::Geom_traits::FT
#endif
seperation_distance(const PointRange1& r1, const PointRange2& r2,
                    const NamedParameters1 np1 = parameters::default_values(),
                    const NamedParameters2 np2 = parameters::default_values());



// TODO: add OBB dedicated code in OBB package

}} // CGAL::Convex_hull_3 namespace

#endif // CGAL_CONVEX_HULL_3_PREDICATES_H
