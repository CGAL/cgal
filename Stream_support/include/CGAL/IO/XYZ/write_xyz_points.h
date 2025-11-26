// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_POINT_SET_PROCESSING_WRITE_XYZ_POINTS_H
#define CGAL_POINT_SET_PROCESSING_WRITE_XYZ_POINTS_H

#include <CGAL/IO/helpers.h>

#include <CGAL/IO/XYZ.h>
#include <CGAL/property_map.h>
#include <CGAL/assertions.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Iterator_range.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iostream>
#include <fstream>
#include <iterator>
#include <type_traits>

namespace CGAL {
namespace Point_set_processing_3 {
namespace internal {

template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_XYZ_PSP(std::ostream& os,
                   const PointRange& points,
                   const CGAL_NP_CLASS& np = CGAL::parameters::default_values())
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  // basic geometric types
  typedef Point_set_processing_3_np_helper<PointRange, CGAL_NP_CLASS> NP_helper;
  typedef typename NP_helper::Const_point_map PointMap;
  typedef typename NP_helper::Normal_map NormalMap;

  const bool has_normals = NP_helper::has_normal_map(points, np);

  PointMap point_map = NP_helper::get_const_point_map(points, np);
  NormalMap normal_map = NP_helper::get_normal_map(points, np);

  CGAL_precondition(points.begin() != points.end());

  if(!os)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  set_stream_precision_from_NP(os, np);

  // Write positions + normals
  for(typename PointRange::const_iterator it = points.begin(); it != points.end(); it++)
  {
    os << get(point_map, *it);
    if(has_normals)
      os << " " << get(normal_map, *it);
    os << "\n";
  }

  os << std::flush;

  return !os.fail();
}

} // namespace internal
} // Point_set_processing_3

namespace IO {

// documented in ../XYZ.h
template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
bool write_XYZ(std::ostream& os,
               const PointRange& points,
               const CGAL_NP_CLASS& np,
               std::enable_if_t<internal::is_Range<PointRange>::value>*
               )
{
  return Point_set_processing_3::internal::write_XYZ_PSP(os, points, np);
}

// documented in ../XYZ.h
template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
bool write_XYZ(const std::string& filename,
               const PointRange& points,
               const CGAL_NP_CLASS& np,
               std::enable_if_t<internal::is_Range<PointRange>::value>*
               )
{
  std::ofstream os(filename);
  return write_XYZ(os, points, np);
}

} // namespace IO

} // namespace CGAL

#endif // CGAL_POINT_SET_PROCESSING_WRITE_XYZ_POINTS_H
