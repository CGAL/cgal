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
// Author(s): Shepard Liu       <shepard0liu@gmail.com>

#ifndef ARR_VIEWER_H
#define ARR_VIEWER_H

#include <algorithm>
#include <array>
#include <cstdlib>
#include <type_traits>

#include <QWidget>

#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Qt/camera.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Basic_viewer.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Qt/camera.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Buffer_for_vao.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Draw_aos/type_utils.h>
#include <CGAL/Draw_aos/Arr_render_context.h>
#include <CGAL/Draw_aos/Arr_bounded_renderer.h>
#include <CGAL/Draw_aos/Arr_coordinate_converter.h>
#include <CGAL/Draw_aos/Arr_face_point_generator.h>

namespace CGAL {
namespace draw_aos {

/*! \brief Viewport helper functions
 *
 * \tparam Arrangement
 */
template <typename Arrangement, typename = void>
class Arr_viewport_helpers;

// Specialization for planar arrangements
template <typename Arrangement>
class Arr_viewport_helpers<Arrangement,
                           std::enable_if_t<! is_or_derived_from_curved_surf_traits_v
                                            <typename Arrangement::Geometry_traits_2>>> {
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Camera = qglviewer::Camera;
  using Point = typename Approx_traits::Point;
  using Local_point = Buffer_for_vao::Local_point;

protected:
  Arr_viewport_helpers(const Arrangement& arr) : m_arr(arr) {}

  /*! \brief Computes a subpixel-level approximation error based on the bounding box and viewport width.
   *
   * \param bbox
   * \param viewport_width width of the viewport in pixels
   * \return double
   */
  double approximation_error(const Bbox_2& bbox, int viewport_width) const
  { return bbox.x_span() / viewport_width; }

  /*! \brief Computes a parameter space bounding box that contains everything in the arrangement with some margin.
   *
   * \note For arrangement induced by unbounded curves, the bounding box only fits all vertices.
   * \return Bbox_2
   */
  Bbox_2 arr_bbox() const {
    const auto& traits = *m_arr.geometry_traits();
    Bbox_2 bbox;
    // Computes a rough bounding box from the vertices.
    for (const auto& vh : m_arr.vertex_handles())
      bbox += traits.approximate_2_object()(vh->point()).bbox();

    double approx_error = approximation_error(bbox, 100);
    // Computes a more precise bounding box from the halfedges.
    auto approx = traits.approximate_2_object();
    for (const auto& he : m_arr.halfedge_handles()) {
      approx(he->curve(), approx_error,
             boost::make_function_output_iterator([&bbox](Approx_point pt) { bbox += pt.bbox(); }));
    }
    // Place margin around the bbox.
    double dx = bbox.x_span() * 0.1;
    double dy = bbox.y_span() * 0.1;
    bbox = Bbox_2(bbox.xmin() - dx, bbox.ymin() - dy, bbox.xmax() + dx, bbox.ymax() + dy);
    // Make sure the bbox is not degenerate.
    if (bbox.x_span() == 0) bbox += Bbox_2(bbox.xmin() - 1, bbox.ymin(), bbox.xmax() + 1, bbox.ymax());
    if (bbox.y_span() == 0) bbox += Bbox_2(bbox.xmin(), bbox.ymin() - 1, bbox.xmax(), bbox.ymax() + 1);
    return bbox;
  }

  /*! \brief Fits the camera to bbox.
   *
   * \param bbox
   * \param camera
   */
  void fit_camera(const Bbox_2& bbox, Camera& cam) const {
    using Vec = qglviewer::Vec;
    cam.fitBoundingBox(Vec(bbox.xmin(), bbox.ymin(), 0.0), Vec(bbox.xmax(), bbox.ymax(), 0.0));
  }

  /*! \brief Computes parameter space axis aligned bounding box from camera parameters.
   *
   * \param cam
   * \return Bbox_2
   */
  Bbox_2 screen_to_world(const Camera& cam) const {
    QMatrix4x4 mvp;
    cam.getModelViewProjectionMatrix(mvp.data());
    QMatrix4x4 inverse_mvp = mvp.inverted();
    // Define 4 corners of the near plane in NDC (-1 to 1 in x and y)
    std::array<QVector4D, 4> clip_space_corners{QVector4D(-1.0, -1.0, 0.0, 1.0), QVector4D(-1.0, 1.0, 0.0, 1.0),
                                                QVector4D(1.0, -1.0, 0.0, 1.0), QVector4D(1.0, 1.0, 0.0, 1.0)};
    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();
    for (const QVector4D& corner : clip_space_corners) {
      QVector4D world = inverse_mvp * corner;
      if (world.w() != 0.0) world /= world.w();
      double x = world.x();
      double y = world.y();
      xmin = std::min(xmin, x);
      xmax = std::max(xmax, x);
      ymin = std::min(ymin, y);
      ymax = std::max(ymax, y);
    }
    return Bbox_2(xmin, ymin, xmax, ymax);
  }

  /*! \brief Converts a parameter space point to a local point of the buffer object.
   *
   * \param pt
   * \return Local_point
   */
  Local_point to_local_point(Point pt) const { return Local_point(pt.x(), pt.y(), 0.0); }

private:
  const Arrangement& m_arr;
};

// Spherical arrangement specialization
template <typename Arrangement>
class Arr_viewport_helpers<Arrangement,
                           std::enable_if_t<is_or_derived_from_agas_v<typename Arrangement::Geometry_traits_2>>> {
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Camera = qglviewer::Camera;
  using Point = typename Approx_traits::Point;
  using Local_point = Buffer_for_vao::Local_point;

protected:
  Arr_viewport_helpers(const Arrangement& arr) : m_arr(arr) {}

  Bbox_2 arr_bbox() const { return Bbox_2(0, 0, 2 * CGAL_PI, CGAL_PI); }

  Bbox_2 screen_to_world(const Camera& /* cam */) const { return Bbox_2(0, 0, 2 * CGAL_PI, CGAL_PI); }

  void fit_camera(const Bbox_2&, Camera& cam) {
    using Vec = qglviewer::Vec;
    cam.setSceneCenter(Vec(0, 0, 0));
    cam.fitSphere(Vec(0, 0, 0), 1.1); // slightly larger than the unit sphere
  }

  double approximation_error(const Bbox_2& bbox, int viewport_width) const {
    // If crossing hemisphere
    if (bbox.x_span() >= CGAL_PI) return 1.0 / viewport_width;
    // Otherwise we evalute the error bound with respect to the longest longitude arc
    double theta =
      std::abs(bbox.ymin() - CGAL_PI / 2.0) < std::abs(bbox.ymax() - CGAL_PI / 2.0) ? bbox.ymin() : bbox.ymax();
    return bbox.x_span() * std::sin(theta) / viewport_width;
  }

  Buffer_for_vao::Local_point to_local_point(Point pt) const {
    auto approx_pt = Arr_coordinate_converter(*m_arr.geometry_traits()).to_cartesian(pt);
    return Buffer_for_vao::Local_point(approx_pt.dx(), approx_pt.dy(), approx_pt.dz());
  }

private:
  const Arrangement& m_arr;
};

/*! Viewer for visualizing arrangements on surface.
 *
 * \tparam Arrangement
 * \tparam GSOptions
 */
template <typename Arrangement, typename GSOptions>
class Arr_viewer : public Qt::Basic_viewer, Arr_viewport_helpers<Arrangement> {
  using Basic_viewer = Qt::Basic_viewer;
  using Helpers = Arr_viewport_helpers<Arrangement>;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Point = typename Approx_traits::Point;
  using Point_generator = Arr_face_point_generator<Arrangement>;
  using Faces_point_map = typename Point_generator::Face_points_map;

  struct Render_params {
    bool operator==(const Render_params& other) const
    { return bbox == other.bbox && approx_error == other.approx_error; }

    Bbox_2 bbox;
    double approx_error{0};
  };

  constexpr static bool Is_on_curved_surface = is_or_derived_from_curved_surf_traits_v<Geom_traits>;

private:
  static bool contains(const Bbox_2& bbox, const Point& pt)
  { return bbox.xmin() <= pt.x() && pt.x() <= bbox.xmax() && bbox.ymin() <= pt.y() && pt.y() <= bbox.ymax(); }

  int viewport_width() const {
    std::array<GLint, 4> viewport;
    this->camera_->getViewport(viewport.data());
    return viewport[2];
  }

  Render_params compute_render_params() {
    Render_params params;
    params.bbox = this->screen_to_world(*this->camera_);
    params.approx_error = this->approximation_error(params.bbox, viewport_width());
    return params;
  }

  void render_arr(const Render_params& params) {
    const Bbox_2& bbox = params.bbox;
    auto face_points = Point_generator()(m_arr, params.approx_error);
    Arr_render_context<Arrangement> ctx(m_arr, params.approx_error, face_points);
    Arr_bounded_renderer<Arrangement> renderer(ctx, bbox);
    auto cache = renderer.render();

    // add faces
    for (const auto& [fh, tf] : cache.faces()) {
      if (! m_gso.draw_face(m_arr, fh)) continue;
      bool colored_face = m_gso.colored_face(m_arr, fh);
      auto color = colored_face ? m_gso.face_color(m_arr, fh) : CGAL::IO::Color();
      for (const auto& tri : tf.triangles) {
        if (colored_face) m_gs.face_begin(color);
        else m_gs.face_begin();
        for (const auto i : tri) m_gs.add_point_in_face(this->to_local_point(tf.points[i]));
        m_gs.face_end();
      }
    }
    // add edges
    for (const auto& [he, polyline] : cache.halfedges()) {
      if (he->direction() == ARR_RIGHT_TO_LEFT || !m_gso.draw_edge(m_arr, he) || polyline.size() < 2) continue;
      bool colored_edge = m_gso.colored_edge(m_arr, he);
      auto color = colored_edge ? m_gso.edge_color(m_arr, he) : CGAL::IO::Color();
      // skip first two if starts with a sep point.
      int start_idx = Approx_traits::is_null(polyline.front()) ? 2 : 0;
      // skip last two if ends with a sep point.
      int end_idx = Approx_traits::is_null(polyline.back()) ? polyline.size() - 2 : polyline.size();
      for (int i = start_idx; i < end_idx - 1; ++i) {
        const auto& src = polyline[i];
        const auto& tgt = polyline[i + 1];
        if (Approx_traits::is_null(src) || Approx_traits::is_null(tgt)) continue;
        if (! contains(bbox, src) || !contains(bbox, tgt)) continue;
        if (colored_edge)
          m_gs.add_segment(this->to_local_point(src), this->to_local_point(tgt), color);
        else
          m_gs.add_segment(this->to_local_point(src), this->to_local_point(tgt));
      }
    }
    // add vertices
    for (const auto& [vh, pt] : cache.vertices()) {
      if (! m_gso.draw_vertex(m_arr, vh) || !contains(bbox, pt)) continue;
      if (m_gso.colored_vertex(m_arr, vh))
        m_gs.add_point(this->to_local_point(pt), m_gso.vertex_color(m_arr, vh));
      else
        m_gs.add_point(this->to_local_point(pt));
    }
  }

  /*! \brief Rerender scene within the given bounding box.
   *
   * \param bbox
   */
  void rerender(const Render_params& params) {
    if (params == m_last_params) return;
    m_last_params = params;
    m_gs.clear();
    render_arr(params);
    Basic_viewer::redraw();
  }

public:
  Arr_viewer(QWidget* parent, const Arrangement& arr, const GSOptions& gso, const char* title, Bbox_2 initial_bbox) :
    Basic_viewer(parent, m_gs, title),
    Helpers(arr),
    m_gso(gso),
    m_arr(arr),
    m_coords(*arr.geometry_traits()) {
    if ((initial_bbox.x_span() == 0) || (initial_bbox.y_span() == 0) || (Is_on_curved_surface))
      m_initial_bbox = this->arr_bbox();
    else
      m_initial_bbox = initial_bbox;
  }

  virtual void draw() override {
    Render_params params = compute_render_params();

#if defined(CGAL_DRAW_AOS_DEBUG)
    if constexpr(! is_or_derived_from_agas_v<Geom_traits>) {
      Bbox_2& bbox = params.bbox;
      double dx = (bbox.xmax() - bbox.xmin()) * 0.1;
      double dy = (bbox.ymax() - bbox.ymin()) * 0.1;
      bbox = Bbox_2(bbox.xmin() + dx, bbox.ymin() + dy, bbox.xmax() - dx, bbox.ymax() - dy);
      std::cout << "Camera changed, recomputing arrangement bounding box: " << bbox << std::endl;
    }
#endif

    rerender(params);

    if (! m_initialized) {
      // The initial render must be done with original camera parameters or the width of edges gets exaggerated.
      // So we fit the camera after initial render.
      this->fit_camera(m_initial_bbox, *this->camera_);
      m_initialized = true;
    }

    Basic_viewer::draw();
  }

  virtual ~Arr_viewer() {}

private:
  Graphics_scene m_gs;
  GSOptions m_gso;
  const Arrangement& m_arr;
  bool m_initialized{false};
  Bbox_2 m_initial_bbox;
  const Arr_coordinate_converter<Geom_traits> m_coords;
  Render_params m_last_params;
};

} // namespace draw_aos
} // namespace CGAL

#endif
