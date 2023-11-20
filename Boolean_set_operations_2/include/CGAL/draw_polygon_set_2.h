// Copyright (c) 1997
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
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_POLYGON_SET_2_H
#define CGAL_DRAW_POLYGON_SET_2_H

#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>

#ifdef DOXYGEN_RUNNING
namespace CGAL {

/*!
 * \ingroup PkgDrawPolygonSet2
 *
 * opens a new window and draws `aps`, an instance of the `CGAL::Polygon_set_2`
 * class. A call to this function is blocking, that is the program continues as
 * soon as the user closes the window. This function requires `CGAL_Qt5`, and is
 * only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.  Linking with
 * the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt5` and add
 * the definition `CGAL_USE_BASIC_VIEWER`.
 * \tparam PS an instance of the `CGAL::Polygon_set_2` class.
 * \param aps the polygon set to draw.
 */
template<class PS>
void draw(const PS& aps);

} /* namespace CGAL */
#endif

#include <CGAL/Polygon_set_2.h>

namespace CGAL {

namespace draw_function_for_boolean_set_2 {

template <typename PS2, class GSOptions>
void compute_loop(const typename PS2::Polygon_2& p, bool hole,
                  CGAL::Graphics_scene& gs,
                  const GSOptions& gso)
{
  if (hole)
  { gs.add_point_in_face(p.vertex(p.size()-1)); }

  auto prev = p.vertices_begin();
  auto it = prev;
  gs.add_point(*it);
  gs.add_point_in_face(*it);
  for (++it; it != p.vertices_end(); ++it)
  {
    gs.add_point(*it);           // add vertex
    gs.add_segment(*prev, *it);  // add segment with previous point
    gs.add_point_in_face(*it);   // add point in face
    prev = it;
  }

  // Add the last segment between the last point and the first one
  gs.add_segment(*prev, *(p.vertices_begin()));
}

/// Compute the elements of a polygon with holes.
template <typename PWH, class GSOptions>
void compute_elements(const PWH& pwh,
                      CGAL::Graphics_scene& gs,
                      const GSOptions& gso)
{
  if (!gso.draw_unbounded() && pwh.outer_boundary().is_empty()) return;

  CGAL::IO::Color c(75,160,255);
  gs.face_begin(c);

  const typename PWH::Point_2* point_in_face;
  if (pwh.outer_boundary().is_empty())
  {
    typename PWH::Polygon_2 pgn;
    pgn.push_back(Pnt(-gso.width(), -gso.height()));
    pgn.push_back(Pnt(gso.width(), -gso.height()));
    pgn.push_back(Pnt(gso.width(), gso.height()));
    pgn.push_back(Pnt(-gso.width(), gso.height()));
    compute_loop(pgn, false, gs);
    point_in_face = &(pgn.vertex(pgn.size()-1));
  }
  else
  {
    const auto& outer_boundary = pwh.outer_boundary();
    compute_loop(outer_boundary, false, gs);
    point_in_face = &(outer_boundary.vertex(outer_boundary.size()-1));
  }

  for (auto it = pwh.holes_begin(); it != pwh.holes_end(); ++it)
  {
    compute_loop(*it, true, gs);
    gs.add_point_in_face(*point_in_face);
  }

  gs.face_end();
}

} // End namespace draw_function_for_boolean_set_2

#ifdef CGAL_USE_BASIC_VIEWER

template <typename PolygonSet_2>
class Polygon_set_2_basic_viewer_qt : public Basic_viewer
{
  using Base = Basic_viewer;
  using Ps = PolygonSet_2;
  using Pwh = typename Ps::Polygon_with_holes_2;
  using Pgn = typename Ps::Polygon_2;
  using Pnt = typename Pgn::Point_2;

public:
  Polygon_set_2_basic_viewer_qt(QWidget* parent, const Ps& ps,
                                const char* title = "Basic Polygon_set_2 Viewer",
                                bool draw_unbounded = false,
                                bool draw_vertices = false) :
    Base(parent, graphics_scene, title, draw_vertices),
    m_ps(ps),
    m_draw_unbounded(draw_unbounded)
  {
    if (ps.is_empty()) return;

    // mimic the computation of Camera::pixelGLRatio()
    auto bbox = bounding_box();
    CGAL::qglviewer::Vec minv(bbox.xmin(), bbox.ymin(), 0);
    CGAL::qglviewer::Vec maxv(bbox.xmax(), bbox.ymax(), 0);
    auto diameter = (maxv - minv).norm();
    m_pixel_ratio = diameter / m_height;
  }

  /*! Intercept the resizing of the window.
   */
  virtual void resizeGL(int width, int height) {
    CGAL::QGLViewer::resizeGL(width, height);
    m_width = width;
    m_height = height;
    CGAL::qglviewer::Vec p;
    auto ratio = this->camera()->pixelGLRatio(p);
    if (ratio != m_pixel_ratio) {
      m_pixel_ratio = ratio;
      add_elements();
    }
  }

  /*! Obtain the pixel ratio.
   */
  double pixel_ratio() const { return m_pixel_ratio; }

  /*! Compute the elements of a polygon set.
   */
  virtual void add_elements()
  {
    graphics_scene.clear();
    std::vector<Pwh> pwhs;
    m_ps.polygons_with_holes(std::back_inserter(pwhs));
    for (const auto& pwh : pwhs)
    { compute_elements(pwh, graphics_scene, graphics_scene_options); }
  }

  /*! Compute the bounding box.
   */
  CGAL::Bbox_2 bounding_box() {
    Bbox_2 bbox;
    std::vector<Pwh> pwhs;
    m_ps.polygons_with_holes(std::back_inserter(pwhs));
    for (const auto& pwh : pwhs) {
      if (! pwh.is_unbounded()) {
        bbox += pwh.outer_boundary().bbox();
        continue;
      }
      for (auto it = pwh.holes_begin(); it != pwh.holes_end(); ++it)
        bbox += it->bbox();
    }
    return bbox;
  }

private:
  //! The window width in pixels.
  int m_width = CGAL_BASIC_VIEWER_INIT_SIZE_X;

  //! The window height in pixels.
  int m_height = CGAL_BASIC_VIEWER_INIT_SIZE_Y;

  //! The ratio between pixel and opengl units (in world coordinate system).
  double m_pixel_ratio = 1;

  //! The polygon set to draw.
  const Ps& m_ps;

  //! Indicates whether to draw unbounded polygons with holes.
  bool m_draw_unbounded = false;

  Graphics_scene graphics_scene;
  Graphics_scene_options<Ps> graphics_scene_options;
};

// Specialization of draw function.
template<class T, class C, class D>
void draw(const CGAL::Polygon_set_2<T, C, D>& ps,
          const char* title = "Polygon_set_2 Basic Viewer",
          bool draw_vertices = false,
          bool draw_unbounded = false)
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite = true;
#else
  bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (! cgal_test_suite)
  {
    using Ps = CGAL::Polygon_set_2<T, C, D>;
    using Viewer = Polygon_set_2_basic_viewer_qt<Ps>;
    CGAL::Qt::init_ogl_context(4,3);
    int argc = 1;
    const char* argv[2] = {"t2_viewer", nullptr};
    QApplication app(argc, const_cast<char**>(argv));
    Viewer basic_viewer(app.activeWindow(), ps,
                        title, draw_unbounded, draw_vertices);
    basic_viewer.add_elements();
    basic_viewer.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POLYGON_SET_2_H
