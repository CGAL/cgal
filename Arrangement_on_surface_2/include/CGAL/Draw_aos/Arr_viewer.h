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

#ifndef ARR_VIEWER_H
#define ARR_VIEWER_H

#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <type_traits>

#include <boost/iterator/function_output_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <QtWidgets/QApplication>
#include <QtWidgets/QWidget>
#include <QtOpenGLWidgets/QtOpenGLWidgets>
#include <QtGui/QOpenGLFunctions>
#include <QtGui/QMouseEvent>
#include <QtGui/QKeyEvent>

#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Qt/camera.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Basic_viewer.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Qt/camera.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Buffer_for_vao.h>
#include <CGAL/Draw_aos/Arr_portals.h>
#include <CGAL/Draw_aos/type_utils.h>
#include <CGAL/Draw_aos/Arr_bounded_renderer.h>
#include <CGAL/Draw_aos/Arr_projector.h>
#include <CGAL/Draw_aos/Arr_render_context.h>

namespace CGAL {
namespace draw_aos {

/**
 * @brief Viewer for visualizing arrangements on surface.
 *
 * @tparam Arrangement should be parameterized by a geometry traits that models Approximate_2, or it silently fails.
 * @tparam GSOptions Graphics scene options
 */
template <typename Arrangement, typename GSOptions, typename = void>
class Arr_viewer;

template <typename Arrangement, typename GSOptions>
class Arr_viewer<Arrangement,
                 GSOptions,
                 std::enable_if_t<!has_approximate_traits_v<typename Arrangement::Geometry_traits_2>>>
    : public Qt::Basic_viewer
{
public:
  Arr_viewer(QWidget* parent, const Arrangement& arr, GSOptions options, const char* title = "Arrangement Viewer")
      : Qt::Basic_viewer(parent, Graphics_scene(), title) {
    std::cerr << "Arr_viewer requires geometry traits that supports approximation of Point_2 and "
                 "X_monotone_curve_2."
              << std::endl;
    exit(1);
  }

  virtual ~Arr_viewer() = default;
};

template <typename Arrangement, typename GSOptions>
class Arr_viewer<Arrangement,
                 GSOptions,
                 std::enable_if_t<has_approximate_traits_v<typename Arrangement::Geometry_traits_2>>>
    : public Qt::Basic_viewer
{
  using Basic_viewer = Qt::Basic_viewer;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Feature_portal_map = typename Arr_portals<Arrangement>::Feature_portals_map;
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Point_geom = typename Approx_traits::Point_geom;

  struct Render_params
  {
    bool operator==(const Render_params& other) const {
      return bbox == other.bbox && approx_error == other.approx_error;
    }
    Bbox_2 bbox;
    double approx_error{0};
  };

private:
  static bool contains(const Bbox_2& bbox, const Point_geom& pt) {
    return bbox.xmin() <= pt.x() && pt.x() <= bbox.xmax() && bbox.ymin() <= pt.y() && pt.y() <= bbox.ymax();
  }

  /**
   * @brief Computes the initial bounding box of the arrangement.
   *
   * @return Bbox_2
   */
  Bbox_2 initial_bbox() const { return initial_bbox_impl<Geom_traits>(); }

  template <typename T, std::enable_if_t<!is_or_derived_from_agas_v<T>, int> = 0>
  Bbox_2 initial_bbox_impl() const {
    const auto& traits = *m_arr.geometry_traits();
    Bbox_2 bbox;
    // Computes a rough bounding box from the vertices.
    for(const auto& vh : m_arr.vertex_handles()) {
      bbox += m_proj.project(traits.approximate_2_object()(vh->point())).bbox();
    }
    double approx_error = get_approx_error(bbox);
    // Computes a more precise bounding box from the halfedges.
    for(const auto& he : m_arr.halfedge_handles()) {
      traits.approximate_2_object()(he->curve(), approx_error,
                                    boost::make_function_output_iterator(
                                        [this, &bbox](Approx_point pt) { bbox += this->m_proj.project(pt).bbox(); }));
    }
    // Place margin around the bbox.
    double dx = bbox.x_span() * 0.1;
    double dy = bbox.y_span() * 0.1;
    bbox = Bbox_2(bbox.xmin() - dx, bbox.ymin() - dy, bbox.xmax() + dx, bbox.ymax() + dy);
    // Make sure the bbox is not degenerate.
    if(bbox.x_span() == 0) bbox += Bbox_2(bbox.xmin() - 1, bbox.ymin(), bbox.xmax() + 1, bbox.ymax());
    if(bbox.y_span() == 0) bbox += Bbox_2(bbox.xmin(), bbox.ymin() - 1, bbox.xmax(), bbox.ymax() + 1);
    return bbox;
  }

  // Specialization for geodesic arc on sphere traits
  template <typename T, std::enable_if_t<is_or_derived_from_agas_v<T>, int> = 0>
  Bbox_2 initial_bbox_impl() const {
    return Bbox_2(0, 0, 2 * CGAL_PI, CGAL_PI);
  }

  /**
   * @brief Computes the corresponding world bounding box of viewport.
   * @return Bbox_2
   */
  Bbox_2 view_bbox_from_camera() const { return view_bbox_from_camera_impl<Geom_traits>(); }

  template <typename T, std::enable_if_t<!is_or_derived_from_agas_v<T>, int> = 0>
  Bbox_2 view_bbox_from_camera_impl() const {
    QMatrix4x4 mvp;
    camera_->getModelViewProjectionMatrix(mvp.data());
    QMatrix4x4 inverse_mvp = mvp.inverted();
    // Define 4 corners of the near plane in NDC (-1 to 1 in x and y)
    std::array<QVector4D, 4> clip_space_corners{QVector4D(-1.0, -1.0, 0.0, 1.0), QVector4D(-1.0, 1.0, 0.0, 1.0),
                                                QVector4D(1.0, -1.0, 0.0, 1.0), QVector4D(1.0, 1.0, 0.0, 1.0)};
    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();
    for(const QVector4D& corner : clip_space_corners) {
      QVector4D world = inverse_mvp * corner;
      if(world.w() != 0.0) world /= world.w();
      double x = world.x();
      double y = world.y();
      xmin = std::min(xmin, x);
      xmax = std::max(xmax, x);
      ymin = std::min(ymin, y);
      ymax = std::max(ymax, y);
    }
    return Bbox_2(xmin, ymin, xmax, ymax);
  }

  // Specialization for geodesic arc on sphere traits
  template <typename T, std::enable_if_t<is_or_derived_from_agas_v<T>, int> = 0>
  Bbox_2 view_bbox_from_camera_impl() const {
    // We fix the bounding box for spherical arrangements to cover the whole sphere.
    return Bbox_2(0, 0, 2 * CGAL_PI, CGAL_PI);
  }

  /**
   * @brief Converts a Point_geom point to euclidean point that can be handled by Buffer_for_vao.
   *
   * @param pt
   * @return Buffer_for_vao::Local_point
   */
  Buffer_for_vao::Local_point transform_point(const Point_geom& proj_pt) const {
    return transform_point_impl<Geom_traits>(proj_pt);
  }

  template <typename T, std::enable_if_t<!is_or_derived_from_agas_v<T>, int> = 0>
  Buffer_for_vao::Local_point transform_point_impl(Point_geom proj_pt) const {
    return Buffer_for_vao::Local_point(proj_pt.x(), proj_pt.y(), 0.0);
  }

  // Specialization for geodesic arc on sphere traits
  template <typename T, std::enable_if_t<is_or_derived_from_agas_v<T>, int> = 0>
  Buffer_for_vao::Local_point transform_point_impl(Point_geom proj_pt) const {
    auto approx_pt = m_proj.unproject(proj_pt);
    return Buffer_for_vao::Local_point(approx_pt.dx(), approx_pt.dy(), approx_pt.dz());
  }

  /**
   * @brief Fits the camera to the given bounding box.
   *
   * Works only for orthogonal camera directed onto a 2D plane by now.
   * @param bbox
   */
  void fit_camera(const Bbox_2& bbox) { fit_camera_impl<Geom_traits>(bbox); }

  template <typename T, std::enable_if_t<!is_or_derived_from_agas_v<T>, int> = 0>
  void fit_camera_impl(const Bbox_2& bbox) {
    // CGAL_assertion(camera_->type() == qglviewer::Camera::ORTHOGRAPHIC);
    using Vec = qglviewer::Vec;
    camera_->fitBoundingBox(Vec(bbox.xmin(), bbox.ymin(), 0.0), Vec(bbox.xmax(), bbox.ymax(), 0.0));
  }

  // Specialization for geodesic arc on sphere traits
  template <typename T, std::enable_if_t<is_or_derived_from_agas_v<T>, int> = 0>
  void fit_camera_impl(const Bbox_2&) {
    using Vec = qglviewer::Vec;
    camera_->setSceneCenter(Vec(0, 0, 0));
    camera_->fitSphere(Vec(0, 0, 0), 1.1); // slightly larger than the sphere
  }

  /**
   * @brief Computes a proper approximation error bound.
   *
   * @param bbox 2D parameter space bounding box
   * @return double
   */
  double get_approx_error(const Bbox_2& bbox) const { return get_approx_error_impl<Geom_traits>(bbox); }

  template <typename T, std::enable_if_t<!is_or_derived_from_agas_v<T>, int> = 0>
  double get_approx_error_impl(const Bbox_2& bbox) const {
    std::array<GLint, 4> viewport;
    camera_->getViewport(viewport.data());
    double viewport_width = static_cast<double>(viewport[2]);
    return bbox.x_span() / viewport_width;
  }

  // Specialization for geodesic arc on sphere traits
  template <typename T, std::enable_if_t<is_or_derived_from_agas_v<T>, int> = 0>
  double get_approx_error_impl(const Bbox_2& bbox) const {
    std::array<GLint, 4> viewport;
    camera_->getViewport(viewport.data());
    double viewport_width = static_cast<double>(viewport[2]);
    // If crossing hemisphere
    if(bbox.x_span() >= CGAL_PI) return 1.0 / viewport_width;
    // Otherwise we evalute the error bound with respect to the longest longitude arc
    double theta =
        std::abs(bbox.ymin() - CGAL_PI / 2.0) < std::abs(bbox.ymax() - CGAL_PI / 2.0) ? bbox.ymin() : bbox.ymax();
    return bbox.x_span() * std::sin(theta) / viewport_width;
  }

  /**
   * @brief Computes the render parameters from camera states.
   */
  Render_params get_arr_render_params() {
    Render_params params;
    params.bbox = view_bbox_from_camera();
    params.approx_error = get_approx_error(params.bbox);
    return params;
  }

  /**
   * @brief Render arrangement within the given bounding box.
   *
   * @param bbox
   */
  void render_arr(const Render_params& params) {
    const Bbox_2& bbox = params.bbox;
    Arr_render_context<Arrangement> ctx(m_arr, m_portals.get(), params.approx_error);
    Arr_bounded_renderer<Arrangement> renderer(ctx, bbox);
    auto cache = renderer.render();
    // add faces
    for(const auto& [fh, face_tris] : cache.faces()) {
      const auto& points = face_tris.points;
      const auto& tris = face_tris.triangles;
      bool draw_face = m_gso.colored_face(m_arr, fh);
      for(const auto& t : tris) {
        if(draw_face)
          m_gs.face_begin(m_gso.face_color(m_arr, fh));
        else
          m_gs.face_begin();
        for(const auto idx : t) m_gs.add_point_in_face(transform_point(points[idx]));
        m_gs.face_end();
      }
    }
    // add edges
    for(const auto& [he, polyline] : cache.halfedges()) {
      if(polyline.size() < 2) continue;
      bool draw_colored_edge = m_gso.colored_edge(m_arr, he);
      auto color = draw_colored_edge ? m_gso.edge_color(m_arr, he) : CGAL::IO::Color();
      for(size_t i = 0; i < polyline.size() - 1; ++i) {
        const auto& cur_pt = polyline[i];
        const auto& next_pt = polyline[i + 1];
        auto mid_pt = CGAL::midpoint(cur_pt, next_pt);
        if(!contains(bbox, mid_pt)) continue;
        if(draw_colored_edge)
          m_gs.add_segment(transform_point(cur_pt), transform_point(next_pt), color);
        else
          m_gs.add_segment(transform_point(cur_pt), transform_point(next_pt));
      }
    }
    // add vertices
    for(const auto& [vh, pt] : cache.vertices()) {
      if(!contains(bbox, pt)) continue;
      if(m_gso.colored_vertex(m_arr, vh))
        m_gs.add_point(transform_point(pt), m_gso.vertex_color(m_arr, vh));
      else
        m_gs.add_point(transform_point(pt));
    }
  }

  /**
   * @brief Rerender scene within the given bounding box.

   * @param bbox
   */
  void rerender(const Render_params& params) {
    if(params == m_last_render_params) return;
    m_last_render_params = params;
    m_gs.clear();
    render_arr(params);
    Basic_viewer::redraw();
  }

public:
  Arr_viewer(QWidget* parent, const Arrangement& arr, GSOptions gso, const char* title = "Arrangement Viewer")
      : Basic_viewer(parent, m_gs, title)
      , m_gso(gso)
      , m_arr(arr)
      , m_proj(*arr.geometry_traits())
      , m_portals(arr) {}

  virtual void draw() override {
    Render_params params = get_arr_render_params();

#if defined(CGAL_DRAW_AOS_DEBUG)
    if constexpr(!is_or_derived_from_agas_v<Geom_traits>) {
      Bbox_2& bbox = params.bbox;
      double dx = (bbox.xmax() - bbox.xmin()) * 0.1;
      double dy = (bbox.ymax() - bbox.ymin()) * 0.1;
      bbox = Bbox_2(bbox.xmin() + dx, bbox.ymin() + dy, bbox.xmax() - dx, bbox.ymax() - dy);
      std::cout << "Camera changed, recomputing arrangement bounding box: " << bbox << std::endl;
    }
#endif

    rerender(params);

    if(!m_initialized) {
      // The initial render must be done with original camera parameters or the width of edges gets exaggerated.
      // So we fit the camera after initial render.
      fit_camera(initial_bbox());
      m_initialized = true;
    }

    Basic_viewer::draw();
  }

  virtual ~Arr_viewer() {}

private:
  Graphics_scene m_gs;
  GSOptions m_gso;
  const Arrangement& m_arr;
  const Arr_portals<Arrangement> m_portals;
  bool m_initialized{false};
  const Arr_projector<Geom_traits> m_proj;
  Render_params m_last_render_params;
};

} // namespace draw_aos
} // namespace CGAL

#endif
