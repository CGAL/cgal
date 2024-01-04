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
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_POLYGON_SET_2_H
#define CGAL_DRAW_POLYGON_SET_2_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef DOXYGEN_RUNNING
namespace CGAL {

/*!
 * \ingroup PkgDrawPolygonSet2
 *
 * opens a new window and draws `aps`, an instance of the `CGAL::Polygon_set_2`
 * class. A call to this function is blocking, that is the program continues as
 * soon as the user closes the window. This function requires `CGAL_Qt6`, and is
 * only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.  Linking with
 * the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add
 * the definition `CGAL_USE_BASIC_VIEWER`.
 * \tparam PS an instance of the `CGAL::Polygon_set_2` class.
 * \param aps the polygon set to draw.
 */
template<class PS>
void draw(const PS& aps);

} /* namespace CGAL */
#endif

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Polygon_set_2.h>

namespace CGAL {

template <typename PolygonSet_2>
class Polygon_set_2_basic_viewer_qt : public Basic_viewer_qt {
  using Base = Basic_viewer_qt;
  using Ps = PolygonSet_2;
  using Pwh = typename Ps::Polygon_with_holes_2;
  using Pgn = typename Ps::Polygon_2;
  using Pnt = typename Pgn::Point_2;

public:
  Polygon_set_2_basic_viewer_qt(QWidget* parent, const Ps& ps,
                                const char* title = "Basic Polygon_set_2 Viewer",
                                bool draw_unbounded = false,
                                bool draw_vertices = false) :
    Base(parent, title, draw_vertices),
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
    auto ratio = camera()->pixelGLRatio(p);
    if (ratio != m_pixel_ratio) {
      m_pixel_ratio = ratio;
      add_elements();
    }
  }

  /*! Compute the elements of a polygon set.
   */
  virtual void add_elements() {
    clear();

    std::vector<Pwh> pwhs;
    m_ps.polygons_with_holes(std::back_inserter(pwhs));
    for (const auto& pwh : pwhs) add_elements(pwh);
  }

  /*! Obtain the pixel ratio.
   */
  double pixel_ratio() const { return m_pixel_ratio; }

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

protected:
  /*! Compute the elements of a polygon with holes.
   */
  void add_elements(const Pwh& pwh) {
    if (! m_draw_unbounded && pwh.outer_boundary().is_empty()) return;

    CGAL::IO::Color c(75,160,255);
    face_begin(c);

    const Pnt* point_in_face;
    if (pwh.outer_boundary().is_empty()) {
      Pgn pgn;
      pgn.push_back(Pnt(-m_width, -m_height));
      pgn.push_back(Pnt(m_width, -m_height));
      pgn.push_back(Pnt(m_width, m_height));
      pgn.push_back(Pnt(-m_width, m_height));
      compute_loop(pgn, false);
      point_in_face = &(pgn.vertex(pgn.size()-1));
    }
    else {
      const auto& outer_boundary = pwh.outer_boundary();
      compute_loop(outer_boundary, false);
      point_in_face = &(outer_boundary.vertex(outer_boundary.size()-1));
    }

    for (auto it = pwh.holes_begin(); it != pwh.holes_end(); ++it) {
      compute_loop(*it, true);
      add_point_in_face(*point_in_face);
    }

    face_end();
  }

  void compute_loop(const Pgn& p, bool hole) {
    if (hole) add_point_in_face(p.vertex(p.size()-1));

    auto prev = p.vertices_begin();
    auto it = prev;
    add_point(*it);
    add_point_in_face(*it);
    for (++it; it != p.vertices_end(); ++it) {
      add_point(*it);           // add vertex
      add_segment(*prev, *it);  // add segment with previous point
      add_point_in_face(*it);   // add point in face
      prev = it;
    }

    // Add the last segment between the last point and the first one
    add_segment(*prev, *(p.vertices_begin()));
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

  if (! cgal_test_suite) {
    using Ps = CGAL::Polygon_set_2<T, C, D>;
    using Viewer = Polygon_set_2_basic_viewer_qt<Ps>;
    CGAL::Qt::init_ogl_context(4,3);
    int argc = 1;
    const char* argv[2] = {"t2_viewer", nullptr};
    QApplication app(argc, const_cast<char**>(argv));
    Viewer mainwindow(app.activeWindow(), ps, title, draw_unbounded, draw_vertices);
    mainwindow.add_elements();
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POLYGON_SET_2_H
