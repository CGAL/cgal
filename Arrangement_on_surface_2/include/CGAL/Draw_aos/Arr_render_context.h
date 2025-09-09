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
// Author(s): Shepard Li        <shepard0liu@gmail.com>

#ifndef CGAL_DRAW_AOS_ARR_RENDER_CONTEXT_H
#define CGAL_DRAW_AOS_ARR_RENDER_CONTEXT_H
#include <cstdlib>
#include <memory>
#include <atomic>
#include <chrono>

#include <CGAL/Bbox_2.h>
#include <CGAL/Arr_point_location_result.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Draw_aos/Arr_approximation_cache.h>
#include <CGAL/Draw_aos/type_utils.h>
#include <CGAL/Draw_aos/Arr_face_point_generator.h>
#include <CGAL/Draw_aos/Arr_coordinate_converter.h>

#if defined(CGAL_DRAW_AOS_DEBUG)
#include <fstream>
#endif

namespace CGAL {
namespace draw_aos {

/**
 * @brief A cancellable context mixin for asynchronous operations. It also tracks elapsed time for performance
 * profiling.
 *
 * The idea is borrowed from golang with a simple implementation.
 * @see https://pkg.go.dev/context
 */
class Arr_cancellable_context_mixin
{
  using Clock = std::chrono::steady_clock;
  using Duration = Clock::duration;
  using Time_point = std::chrono::time_point<Clock, Duration>;

protected:
  Arr_cancellable_context_mixin()
      : m_start_time(Clock::now())
      , m_cancelled(std::make_shared<std::atomic<bool>>(false)) {}

public:
  Time_point start_time() const { return m_start_time; }
  Time_point end_time() const { return m_end_time; }
  Duration elapsed_time() const { return Clock::now() - m_start_time; }
  bool is_cancelled() const { return m_cancelled->load(); }

  void cancel() {
    m_cancelled->store(true, std::memory_order_relaxed);
    m_end_time = Clock::now();
  }

private:
  Time_point m_start_time, m_end_time;
  std::shared_ptr<std::atomic<bool>> m_cancelled;
};

/**
 * @brief Boundary context mixin for rendering arrangements within a bounding box.
 * Provides extended functionality for checking point-bbox relations.
 *
 * @tparam GeomTraits the geometry traits class.
 */
template <typename GeomTraits>
class Arr_bounds_context_mixin
{
  using Geom_traits = GeomTraits;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Point = typename Approx_traits::Point;
  using Approx_nt = typename Approx_traits::Approx_nt;

protected:
  Arr_bounds_context_mixin(const Bbox_2& bbox)
      : m_bbox(bbox) {}

public:
  double xmin() const { return m_bbox.xmin(); }
  double xmax() const { return m_bbox.xmax(); }
  double ymin() const { return m_bbox.ymin(); }
  double ymax() const { return m_bbox.ymax(); }
  const Bbox_2& bbox() const { return m_bbox; }

  bool contains_x(Approx_nt x) const { return xmin() <= x && x <= xmax(); }
  bool contains_y(Approx_nt y) const { return ymin() <= y && y <= ymax(); }
  bool contains(Point pt) const { return contains_x(pt.x()) && contains_y(pt.y()); }

  Point top_left() const { return Point(xmin(), ymax()); }
  Point top_right() const { return Point(xmax(), ymax()); }
  Point bottom_left() const { return Point(xmin(), ymin()); }
  Point bottom_right() const { return Point(xmax(), ymin()); }

  bool is_on_left(Point pt) const { return pt.x() == xmin() && contains_y(pt.y()); }
  bool is_on_right(Point pt) const { return pt.x() == xmax() && contains_y(pt.y()); }
  bool is_on_bottom(Point pt) const { return pt.y() == ymin() && contains_x(pt.x()); }
  bool is_on_top(Point pt) const { return pt.y() == ymax() && contains_x(pt.x()); }
  bool is_on_boundary(Point pt) const { return is_on_left(pt) || is_on_right(pt) || is_on_bottom(pt) || is_on_top(pt); }

private:
  const Bbox_2 m_bbox;
};

template <typename GeomTraits>
using Arr_parameterization_context_mixin = Arr_coordinate_converter<GeomTraits>;

template <typename Arrangement>
class Arr_render_context : public Arr_cancellable_context_mixin,
                           public Arr_parameterization_context_mixin<typename Arrangement::Geometry_traits_2>
{
  using Cancellable_context_mixin = Arr_cancellable_context_mixin;
  using Param_context_mixin = Arr_parameterization_context_mixin<typename Arrangement::Geometry_traits_2>;
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Face_points_map = typename Arr_face_point_generator<Arrangement>::Face_points_map;

public:
  Arr_render_context(const Arrangement& arr, double approx_error, Face_points_map& face_points)
      : Cancellable_context_mixin()
      , Param_context_mixin(*arr.geometry_traits())
      , m_arr(arr)
      , m_traits(*arr.geometry_traits())
      , m_approx_error(approx_error)
      , m_face_points(face_points) {
#if defined(CGAL_DRAW_AOS_DEBUG) && defined(CGAL_DRAW_AOS_TRIANGULATOR_DEBUG_FILE_DIR)
    std::filesystem::path debug_file_dir(CGAL_DRAW_AOS_TRIANGULATOR_DEBUG_FILE_DIR);
    // clear the index file.
    std::filesystem::remove(debug_file_dir / "index.txt");
#endif
  }

public:
  const double m_approx_error;
  const Arrangement& m_arr;
  const Geom_traits& m_traits;
  const Face_points_map& m_face_points;

#if defined(CGAL_DRAW_AOS_DEBUG)
  std::shared_ptr<int> debug_counter = std::make_shared<int>(0);
#endif
};

template <typename Arrangement>
class Arr_bounded_render_context : public Arr_render_context<Arrangement>,
                                   public Arr_bounds_context_mixin<typename Arrangement::Geometry_traits_2>
{
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_point = typename Geom_traits::Approximate_point_2;
  using Render_context = Arr_render_context<Arrangement>;
  using Bounds_context_mixin = Arr_bounds_context_mixin<Geom_traits>;
  using Approx_cache = Arr_approximation_cache<Arrangement>;

public:
  Arr_bounded_render_context(const Render_context& ctx, const Bbox_2& bbox, Approx_cache& cache)
      : Render_context(ctx)
      , Bounds_context_mixin(bbox)
      , m_cache(cache) {}

public:
  Approx_cache& m_cache;
};

} // namespace draw_aos
} // namespace CGAL
#endif
