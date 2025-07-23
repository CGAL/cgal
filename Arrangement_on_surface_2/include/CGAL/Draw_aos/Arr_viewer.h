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
#include <CGAL/Draw_aos/type_utils.h>
#include <CGAL/Draw_aos/Arr_bounded_renderer.h>
#include <CGAL/Draw_aos/Arr_render_context.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Qt/camera.h>

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
  using Point_location = Arr_trapezoid_ric_point_location<Arrangement>;
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_point = typename Arr_approximate_traits<Geom_traits>::Approx_point;

private:
  static bool contains(const Bbox_2& bbox, const Approx_point& pt) {
    return bbox.xmin() <= pt.x() && pt.x() <= bbox.xmax() && bbox.ymin() <= pt.y() && pt.y() <= bbox.ymax();
  }

  // Function to check if the camera's state has changed
  bool is_camera_changed() {
    QMatrix4x4 proj_mat, mv_mat;
    this->camera_->computeProjectionMatrix();
    this->camera_->computeModelViewMatrix();
    this->camera_->getProjectionMatrix(proj_mat.data());
    this->camera_->getModelViewMatrix(mv_mat.data());
    if(proj_mat == m_last_proj_matrix && mv_mat == m_last_modelview_matrix) return false;
    m_last_proj_matrix = proj_mat;
    m_last_modelview_matrix = mv_mat;
    return true;
  }

  void fill_background(Bbox_2 bbox, CGAL::IO::Color color = CGAL::IO::Color(255, 255, 255)) {
    m_gs.face_begin(color);
    m_gs.add_point_in_face(Approx_point(bbox.xmin(), bbox.ymin()));
    m_gs.add_point_in_face(Approx_point(bbox.xmax(), bbox.ymin()));
    m_gs.add_point_in_face(Approx_point(bbox.xmax(), bbox.ymax()));
    m_gs.add_point_in_face(Approx_point(bbox.xmin(), bbox.ymax()));
    m_gs.face_end();
  }

  Bbox_2 initial_bbox() const {
    Bbox_2 bbox;
    // Computes a rough bounding box from the vertices.
    for(const auto& vh : m_arr.vertex_handles()) {
      Approx_point pt = m_arr.geometry_traits()->approximate_2_object()(vh->point());
      bbox += pt.bbox();
    }
    // Make sure the bbox is not degenerate.
    if(bbox.x_span() == 0) bbox += Bbox_2(bbox.xmin() - 1, bbox.ymin(), bbox.xmax() + 1, bbox.ymax());
    if(bbox.y_span() == 0) bbox += Bbox_2(bbox.xmin(), bbox.ymin() - 1, bbox.xmax(), bbox.ymax() + 1);
    return bbox;
  }

  /**
   * @brief Computes the bounding box of the view from orthogonal camera.
   * @return Bbox_2
   */
  Bbox_2 view_bbox_from_camera() const {
    QMatrix4x4 modelview_matrix, projection_matrix;
    camera_->getModelViewMatrix(modelview_matrix.data());
    camera_->getProjectionMatrix(projection_matrix.data());
    QMatrix4x4 inverse_mvp = (projection_matrix * modelview_matrix).inverted();
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

  /**
   * @brief Fits the camera to the given bounding box.
   *
   * Works only for orthogonal camera directed onto a 2D plane by now.
   * @param bbox
   */
  void fit_camera(Bbox_2 bbox) {
    // CGAL_assertion(camera_->type() == qglviewer::Camera::ORTHOGRAPHIC);
    using Vec = qglviewer::Vec;
    camera_->fitBoundingBox(Vec(bbox.xmin(), bbox.ymin(), 0.0), Vec(bbox.xmax(), bbox.ymax(), 0.0));
  }

  double get_approx_error(const Bbox_2& bbox) const {
    std::array<GLint, 4> viewport;
    camera_->getViewport(viewport.data());
    double viewport_width = static_cast<double>(viewport[2]);
    return bbox.x_span() / viewport_width;
  }

  /*!
   */
  void render_arr(const Bbox_2& bbox) {
    Arr_render_context<Arrangement> ctx(m_arr, m_pl, m_feature_portals, get_approx_error(bbox));
    Arr_bounded_renderer<Arrangement> renderer(ctx, bbox);
    const auto& cache = renderer.render();

#if defined(CGAL_DRAW_AOS_DEBUG)
    std::cout << "Rendering arrangement with " << cache.face_cache_size() << " visible faces, "
              << cache.halfedge_cache_size() << " visible halfedges, " << cache.vertex_cache_size()
              << " visible vertices." << std::endl;
#endif

    // add faces
    for (const auto& [fh, face_tris] : cache.face_cache()) {
      const auto& points = face_tris.points;
      const auto& tris = face_tris.triangles;
      bool draw_face = m_gso.colored_face(m_arr, fh);
      for(const auto& t : tris) {
        if(draw_face)
          m_gs.face_begin(m_gso.face_color(m_arr, fh));
        else
          m_gs.face_begin();
        for(const auto idx : t) m_gs.add_point_in_face(points[idx]);
        m_gs.face_end();
      }
    }

    // add edges
    for(const auto& [he, polyline] : cache.halfedge_cache()) {
      if(polyline.size() < 2) continue;
      bool draw_colored_edge = m_gso.colored_edge(m_arr, he);
      auto color = draw_colored_edge ? m_gso.edge_color(m_arr, he) : CGAL::IO::Color();
      for(size_t i = 0; i < polyline.size() - 1; ++i) {
        const auto& cur_pt = polyline[i];
        const auto& next_pt = polyline[i + 1];
        auto mid_pt = CGAL::midpoint(cur_pt, next_pt);
        if(!contains(bbox, mid_pt)) continue;
        if(draw_colored_edge)
          m_gs.add_segment(cur_pt, next_pt, color);
        else
          m_gs.add_segment(cur_pt, next_pt);
      }
    }

    // add vertices
    for(const auto& [vh, pt] : cache.vertex_cache()) {
      if(!contains(bbox, pt)) continue;
      if(m_gso.colored_vertex(m_arr, vh))
        m_gs.add_point(pt, m_gso.vertex_color(m_arr, vh));
      else
        m_gs.add_point(pt);
    }

    // keep scene non-empty to make sure that the Basic_viewer works in 2D mode for planar arrangements.
    if(m_gs.empty()) fill_background(bbox);
  }

  /*!
   */
  void rerender(Bbox_2 bbox) {
    m_gs.clear();
    render_arr(bbox);
    Basic_viewer::redraw();
  }

public:
  Arr_viewer(QWidget* parent, const Arrangement& arr, GSOptions gso, const char* title = "Arrangement Viewer")
      : Basic_viewer(parent, m_gs, title)
      , m_gso(gso)
      , m_arr(arr)
      , m_feature_portals(Arr_portals<Arrangement>().create(arr))
      , m_pl(arr) {}

  virtual void draw() override {
    if (is_camera_changed()) {
      Bbox_2 bbox = view_bbox_from_camera();

#if defined(CGAL_DRAW_AOS_DEBUG)
      double dx = (bbox.xmax() - bbox.xmin()) * 0.1;
      double dy = (bbox.ymax() - bbox.ymin()) * 0.1;
      bbox = Bbox_2(bbox.xmin() + dx, bbox.ymin() + dy, bbox.xmax() - dx, bbox.ymax() - dy);
      std::cout << "Camera changed, recomputing arrangement bounding box: " << bbox << std::endl;
#endif

      rerender(bbox);
    }

    if(!m_initialized) {
      // The initial render must be done with original camera parameters or the width of edges gets exaggerated.
      // So we fit the camera after initial render.
      fit_camera(initial_bbox());
      m_initialized = true;
    }

    Basic_viewer::draw();
  }

  /*!
   */
  virtual ~Arr_viewer() {}

private:
  Graphics_scene m_gs;
  GSOptions m_gso;
  const Arrangement& m_arr;
  Point_location m_pl;
  Feature_portal_map m_feature_portals;
  QMatrix4x4 m_last_proj_matrix;
  QMatrix4x4 m_last_modelview_matrix;
  bool m_initialized{false};
};

} // namespace draw_aos
} // namespace CGAL

#endif
