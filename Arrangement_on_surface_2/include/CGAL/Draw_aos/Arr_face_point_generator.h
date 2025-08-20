// Copyright (c) 2025
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Shepard Liu	 <shepard0liu@gmail.com>

#ifndef CGAL_DRAW_AOS_ARR_FACE_POINT_GENERATOR_H
#define CGAL_DRAW_AOS_ARR_FACE_POINT_GENERATOR_H
#include <utility>
#include <variant>
#include <vector>

#include <boost/iterator/function_output_iterator.hpp>

#include "CGAL/unordered_flat_map.h"
#include "CGAL/Arr_batched_point_location.h"
#include "CGAL/Arr_point_location_result.h"
#include "CGAL/Draw_aos/Arr_coordinate_converter.h"
#include "CGAL/Draw_aos/type_utils.h"

namespace CGAL {
namespace draw_aos {

/*!
 * \brief Generate face interior points.
 *
 * \tparam Arrangement
 */
template <typename Arrangement, typename = void>
class Arr_face_point_generator;

template <typename Arrangement>
class Arr_face_point_generator<
    Arrangement,
    std::enable_if_t<!is_or_derived_from_curved_surf_traits_v<typename Arrangement::Geometry_traits_2>>>
{
  using Point_geom = typename Arr_approximate_traits<typename Arrangement::Geometry_traits_2>::Point;
  using Face_const_handle = typename Arrangement::Face_const_handle;

public:
  using Face_points_map = unordered_flat_map<Face_const_handle, std::vector<Point_geom>>;

  // No-op implementation for non-curved surface arrangements.
  Face_points_map operator()(const Arrangement&, double) { return {}; }
};

template <typename Arrangement>
class Arr_face_point_generator<Arrangement,
                               std::enable_if_t<is_or_derived_from_agas_v<typename Arrangement::Geometry_traits_2>>>
{
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_nt = typename Approx_traits::Approx_nt;
  using Point = typename Approx_traits::Point;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Gt_point = typename Geom_traits::Point_2;
  using Query_result = std::pair<Gt_point, typename Arr_point_location_result<Arrangement>::Type>;

public:
  using Face_points_map = unordered_flat_map<Face_const_handle, std::vector<Point>>;

  Face_points_map operator()(const Arrangement& arr, double error) {
    const Geom_traits& traits = *arr.geometry_traits();

    // Grid sampling in parameter space.
    Approx_nt cell_size = 2.0 * std::acos(1 - error);
    std::vector<Gt_point> points;
    Arr_coordinate_converter<Geom_traits> coords(traits);
    points.reserve(2 * CGAL_PI / cell_size * CGAL_PI / cell_size);
    for(Approx_nt x = 0; x < 2 * CGAL_PI; x += cell_size) {
      for(Approx_nt y = 0; y < CGAL_PI; y += cell_size) {
        auto pt = coords.to_cartesian(Point(x, y));
        points.push_back(traits.construct_point_2_object()(pt.dx(), pt.dy(), pt.dz()));
      }
    }

    unordered_flat_map<Face_const_handle, std::vector<Point>> face_points;
    CGAL::locate(arr, points.begin(), points.end(),
                 boost::make_function_output_iterator([&face_points, &traits, &coords](const Query_result& res) {
                   if(!std::holds_alternative<Face_const_handle>(res.second)) return;
                   Face_const_handle fh = std::get<Face_const_handle>(res.second);
                   auto [it, _] = face_points.try_emplace(fh, std::vector<Point>());
                   it->second.push_back(coords.to_uv(traits.approximate_2_object()(res.first)));
                 }));
    return face_points;
  }
};

} // namespace draw_aos
} // namespace CGAL

#endif