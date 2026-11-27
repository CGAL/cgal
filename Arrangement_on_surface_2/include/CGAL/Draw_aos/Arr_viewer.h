// Copyright (c) 2025 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Shepard Liu       <shepard0liu@gmail.com>

#ifndef CGAL_DRAW_AOS_ARR_VIEWER_H
#define CGAL_DRAW_AOS_ARR_VIEWER_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <algorithm>
#include <array>
#include <limits>
#include <type_traits>

#include <boost/utility/base_from_member.hpp>

#include <QMatrix4x4>
#include <QVector4D>
#include <QWidget>

#include <CGAL/Bbox_2.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/number_type_config.h>
#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Qt/camera.h>
#include <CGAL/Draw_aos/Arr_graphics_scene_builder.h>
#include <CGAL/Draw_aos/type_utils.h>

namespace CGAL {
namespace draw_aos {

/*! \brief Camera helpers, which depend on the parameter space.
 *
 * \tparam Geom_traits
 */
template <typename GeomTraits, typename = void>
struct Arr_camera_helpers;

// Specialization for planar arrangements
template <typename GeomTraits>
struct Arr_camera_helpers<GeomTraits, std::enable_if_t<! is_or_derived_from_curved_surf_traits_v<GeomTraits>>> {
  using Camera = qglviewer::Camera;

  //! fits the camera to bbox.
  static void fit_camera(const Bbox_2& bbox, Camera& cam) {
    using Vec = qglviewer::Vec;
    cam.fitBoundingBox(Vec(bbox.xmin(), bbox.ymin(), 0.0), Vec(bbox.xmax(), bbox.ymax(), 0.0));
  }

  //! computes the parameter space axis-aligned bounding box visible by the camera.
  static Bbox_2 screen_to_world(const Camera& cam) {
    QMatrix4x4 mvp;
    cam.getModelViewProjectionMatrix(mvp.data());
    QMatrix4x4 inverse_mvp = mvp.inverted();
    // Define 4 corners of the near plane in NDC (-1 to 1 in x and y)
    std::array<QVector4D, 4> clip_space_corners{QVector4D(-1.0, -1.0, 0.0, 1.0), QVector4D(-1.0, 1.0, 0.0, 1.0),
                                                QVector4D(1.0, -1.0, 0.0, 1.0), QVector4D(1.0, 1.0, 0.0, 1.0)};
    double xmin = (std::numeric_limits<double>::max)();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = (std::numeric_limits<double>::max)();
    double ymax = std::numeric_limits<double>::lowest();
    for (const QVector4D& corner : clip_space_corners) {
      QVector4D world = inverse_mvp * corner;
      if (world.w() != 0.0) world /= world.w();
      double x = world.x();
      double y = world.y();
      xmin = (std::min)(xmin, x);
      xmax = (std::max)(xmax, x);
      ymin = (std::min)(ymin, y);
      ymax = (std::max)(ymax, y);
    }
    return Bbox_2(xmin, ymin, xmax, ymax);
  }
};

// Specialization for arrangements on the sphere
template <typename GeomTraits>
struct Arr_camera_helpers<GeomTraits, std::enable_if_t<is_or_derived_from_agas_v<GeomTraits>>> {
  using Camera = qglviewer::Camera;

  //! fits the camera to the unit sphere.
  static void fit_camera(const Bbox_2&, Camera& cam) {
    using Vec = qglviewer::Vec;
    cam.setSceneCenter(Vec(0, 0, 0));
    cam.fitSphere(Vec(0, 0, 0), 1.1); // slightly larger than the unit sphere
  }

  //! obtains the entire parameter space; the entire sphere is always rendered.
  static Bbox_2 screen_to_world(const Camera&) { return Bbox_2(0, 0, 2 * CGAL_PI, CGAL_PI); }
};

/*! Viewer for visualizing arrangements on surface.
 *
 * The scene is rebuilt whenever the visible part of the parameter space or the
 * resolution changes, so that the curves are approximated with an error below
 * the pixel size.
 *
 * \tparam Arrangement
 * \tparam GSOptions
 */
template <typename Arrangement, typename GSOptions>
class Arr_viewer :
    // The scene must be constructed before the base viewer, which is given a reference to it.
    private boost::base_from_member<Graphics_scene>,
    public Qt::Basic_viewer {
  using Base_from_member = boost::base_from_member<Graphics_scene>;
  using Basic_viewer = Qt::Basic_viewer;
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Camera_helpers = Arr_camera_helpers<Geom_traits>;
  using Scene_builder = Arr_graphics_scene_builder<Arrangement, GSOptions>;

  struct Render_params {
    bool operator==(const Render_params& other) const
    { return bbox == other.bbox && approx_error == other.approx_error; }

    Bbox_2 bbox;
    double approx_error{0};
  };

  constexpr static bool Is_on_curved_surface = is_or_derived_from_curved_surf_traits_v<Geom_traits>;

  Render_params compute_render_params() const {
    Render_params params;
    params.bbox = Camera_helpers::screen_to_world(*this->camera());
    params.approx_error = m_builder.approximation_error(params.bbox, this->camera()->screenWidth());
    return params;
  }

  //! rebuilds the scene if the render parameters have changed.
  void rerender(const Render_params& params) {
    if (params == m_last_params) return;
    m_last_params = params;
    m_gs.clear();
    m_builder(m_gs, params.bbox, params.approx_error);
    Basic_viewer::redraw();
  }

public:
  /*! constructs.
   *
   * \param initial_bbox the parameter space bounding box shown initially. If
   *        it is degenerate (or the arrangement is on a curved surface), a
   *        bounding box that contains the arrangement is used.
   */
  Arr_viewer(QWidget* parent, const Arrangement& arr, const GSOptions& gso, const char* title, Bbox_2 initial_bbox) :
    Base_from_member(),
    Basic_viewer(parent, Base_from_member::member, title),
    m_gs(Base_from_member::member),
    m_builder(arr, gso) {
    if ((initial_bbox.x_span() <= 0) || (initial_bbox.y_span() <= 0) || (Is_on_curved_surface))
      m_initial_bbox = m_builder.bbox();
    else
      m_initial_bbox = initial_bbox;
  }

  virtual void draw() override {
    Render_params params = compute_render_params();

#if defined(CGAL_DRAW_AOS_DEBUG)
    if constexpr(! Is_on_curved_surface) {
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
      Camera_helpers::fit_camera(m_initial_bbox, *this->camera());
      m_initialized = true;
    }

    Basic_viewer::draw();
  }

private:
  Graphics_scene& m_gs;
  Scene_builder m_builder;
  bool m_initialized{false};
  Bbox_2 m_initial_bbox;
  Render_params m_last_params;
};

} // namespace draw_aos
} // namespace CGAL

#endif
